module system
 implicit none

 integer,parameter :: q=3 !only for 3-state Potts model
 integer :: n,l,step1,step2,bins,step3
 integer, allocatable :: spin(:)
 integer, allocatable :: nbor(:,:)
 integer :: enbor, vect       
 real(8) ::  bp
 real(8), allocatable :: m_dim(:),clusters(:)
 real(8) :: temp,arm=0.d0,asm=0.d0,m_4=0.d0,binder=0.d0,engy=0.d0,engy2=0.d0
 real(8) :: test0val
 real(8),parameter :: pi2=2.d0*3.14159265358979323846
 real(8),parameter :: coff=sqrt(3.d0)/2.d0
 
end module system

program main
 use system
 implicit none
 integer :: i,j,k
  call initran(1)
  open(1,file='read.in',status='old')
  read(1,*)l,temp,step1,step2,bins
  close(1)
  call initial()
  do i=1,step1
    call wolff()
  !  call mcstep()
  enddo
  do i=1,step2        
    do j=1,bins           
      call wolff()
    !  call mcstep()
      call sample  
    enddo                                                             
    call writers
    call checkpoint()
  enddo
  call record()

 deallocate(nbor,spin,clusters)
end program

subroutine initial()
	use system
	implicit none
	integer i,j,x,y,z
  logical :: spin_exist
  real(8) :: norm,rn

  n=l*l
  allocate (clusters(n))
  allocate (nbor(0:n-1,4))
  allocate (spin(0:n-1))

  inquire(file='spin.dat',exist=spin_exist)
  if(.false.) then
      call readin()
    else
      do i=0,n-1   
        spin(i)=int(dble(q)*rn())
      enddo
  endif

  bp=1-exp(-1.d0/temp)
  spin(0)=-1
  nbor=0

  do i=0,n-1
  !2d square lattice with periodic boundary condition
    x=mod(i,l); y=i/l          ! coordinates of site s
    nbor(i,1)=mod(x+1,l)+y*l   ! spin at right neighbor of site s
    nbor(i,2)=x+mod(y+1,l)*l   ! up
    nbor(i,3)=mod(x-1+l,l)+y*l ! left 
    nbor(i,4)=x+mod(y-1+l,l)*l ! down
  end do

end subroutine initial

subroutine writers
  use system
  implicit none
  arm=arm/dble(bins)
  asm=asm/dble(bins)
  m_4=m_4/dble(bins)
  binder=m_4/asm**2
  engy=engy/dble(bins)
  engy2=engy2/dble(bins)

  open(2,file='out.dat',position='append')    
  write(2,1)l,temp,arm,asm,m_4,binder,engy,engy2,test0val,step1,step2,bins
  1 format(i9,'  ',f7.5,7f18.12,3i10)
  close(2)

  arm=0;asm=0;m_4=0;binder=0;
  engy=0;engy2=0;
  test0val=0;
end subroutine writers

subroutine sample
  use system
  implicit none
 	  real(8) :: comb,ans,m,m1,m2,e=0
    integer :: i,j,k
    integer(4) :: n0,n1,n2
    real(8) :: mx,my
    n1=count(spin==1)
    n2=count(spin==2)
    n0=n-n1-n2
    mx=(n0-(n1+n2)/2.d0)/dble(n)
    my=coff*(n1-n2)/dble(n)
    m=mx**2+my**2
    asm=asm+m 
    m_4=m_4+m**2 
    arm=arm+sqrt(m)

    do i=0,n-1
      do j=1,4
        k=nbor(i,j)
        ! if (k.eq.0) cycle
        if (spin(k).eq.spin(i)) e=e-1
      enddo
    enddo
    e=e/(dble(n)*2.d0)
    engy=engy+e
    engy2=engy2+e**2

end subroutine

subroutine mcstep()
  use system
  implicit none
  integer :: i,j,k,s
  real(8) :: rn,ei,ef,ans

  do s=0,n-1
    i=int(rn()*n)
    enbor=0
    vect=int(q*rn())
    if (vect.eq.spin(i)) cycle
    do j=1,4
      k=nbor(i,j)
      ! if (k.eq.0) cycle
      if (spin(k).eq.spin(i)) enbor=enbor-1
      if (spin(k).eq.vect) enbor=enbor+1
    enddo
    if (rn()<=exp(1.d0/temp*enbor)) spin(i)=vect
  enddo
    
endsubroutine

subroutine wolff()
  use system
  implicit none
  integer :: i,j,k,flip,s0,sc,cluster=0
  real(8) :: rn
  i=int(rn()*n)
  s0=spin(i)
  flip=s0+int(rn()*2)+1
  if (flip > q-1) flip=flip-q
  spin(i)=flip
  sc=i
  clusters=0
  do 
    do j=1,4
      k=nbor(sc,j)
      ! if (k.eq.0) cycle
      if ((rn().lt.bp) .and. (spin(k).eq.s0)) then
          spin(k)=flip
          cluster=cluster+1
          clusters(cluster)=k
      endif
    enddo
    if(cluster==0) exit
    sc=clusters(cluster)
    cluster=cluster-1
  enddo
endsubroutine

subroutine record()
  use system
  open(12,file='spin.dat')
  write(12,*)spin
  close(12)
end subroutine record
subroutine readin()
  use system
  open(12,file='spin.dat',status='old')
  read(12,*)spin
  close(12)
end subroutine readin

subroutine checkpoint()
  logical :: stopfile
  inquire(file='stop',exist=stopfile)
  if(stopfile) then
    call record()
    stop
  endif
endsubroutine checkpoint

real(8) function rn()
!-----------------------------------------------------!
! 64-bit linear congruental random number generator   !
! iran64=oran64*2862933555777941757+1013904243        !
!-----------------------------------------------------!
 implicit none

 real(8)    :: dmu64
 integer(8) :: ran64,mul64,add64
 common/bran64/dmu64,ran64,mul64,add64

 ran64=ran64*mul64+add64
 rn=0.5d0+dmu64*dble(ran64)

 end function rn

subroutine initran(w)

implicit none

 integer(8) :: irmax
 integer(4) :: w,nb,b

 real(8)    :: dmu64
 integer(8) :: ran64,mul64,add64
 common/bran64/dmu64,ran64,mul64,add64
      
 irmax=2_8**31
 irmax=2*(irmax**2-1)+1
 mul64=2862933555777941757_8
 add64=1013904243
 dmu64=0.5d0/dble(irmax)

 open(10,file='seed.in',status='old')
 read(10,*)ran64
 close(10)
 if (w.ne.0) then
    open(10,file='seed.in',status='unknown')
    write(10,*)abs((ran64*mul64)/5+5265361)
    close(10)
 endif

 end subroutine initran

