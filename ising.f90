! 需构建邻接表 nbor.txt
module system
 implicit none
 integer,parameter :: o_n=1
 integer,parameter :: cor=20
 integer :: n,step1,step2,bins,step3
 integer, allocatable :: spin(:)
 integer, allocatable :: nbor(:,:)
 integer :: enbor, vect       
 real(8) ::  bp
 real(8), allocatable :: m_dim(:),clusters(:)
 real(8) :: temp,arm,asm,m_4,binder,engy,engy2
 real(8) :: test0val
 real(8),parameter :: pi=dacos(-1.d0)
 
end module system

program main
 use system
 implicit none
 integer :: i,j,k
  call initran(1)
  open(1,file='read.in',status='old')
  read(1,*)n,temp,step1,step2,bins
  close(1)
  call initial()
  do i=1,step1
    call wolff()
  !  call mcstep()
  enddo
  do i=1,step2
    call cleardata         
    do j=1,bins           
      call wolff()
   !   call mcstep()
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

  allocate (clusters(n))
  allocate (nbor(n,cor))
  allocate (spin(0:n))

  inquire(file='spin.dat',exist=spin_exist)
  if(.false.) then
      call readin()
    else
      do i=1,n   
        spin(i)=2*int(2.d0*rn())-1
      enddo
  endif

  bp=1-exp(-2.d0/temp)
  spin(0)=0
  nbor=0
  open(21,file='nbor.txt',status='old')
  do i=1,n
      read(21,*) nbor(i,:)
  end do
  print*,nbor(1,:)
  close(21)

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
  write(2,1)n,temp,arm,asm,m_4,binder,engy,engy2,test0val,step1,step2,bins
  1 format(i9,'  ',f7.5,7f18.12,3i10)
  close(2)

end subroutine writers

subroutine sample
  use system
  implicit none
 	  real(8) :: comb,ans,m,m1,m2,e=0
    integer :: i,j,k
    integer :: x,y,z,x1,y1,z1
    integer :: ans1
    m=sum(spin(:))/dble(n)
    arm=arm+abs(m)
    m=m**2
    asm=asm+m 
    m_4=m_4+m**2 
    do j=1,n
      enbor=0
      do k=1,4
        enbor=enbor+spin(nbor(j,k))
      enddo
      e=e+spin(j)*enbor
    enddo
    e=e/(dble(n)*2.d0)
    engy=engy-e
    engy2=engy2+e**2

end subroutine

subroutine mcstep()
  use system
  implicit none
  integer :: i,j,k,s
  real(8) :: rn,ei,ef,ans

  do s=1,n
    i=int(rn()*n)+1
    enbor=0
    vect=int(2.d0*rn())*2-1
    if (vect.eq.spin(i)) cycle
    do j=1,cor
        enbor=enbor+spin(nbor(i,j))
    enddo
    ei=spin(i)*enbor
    ef=vect*enbor
    if (rn()<=exp((ef-ei)/temp)) spin(i)=vect
  enddo
    
endsubroutine

subroutine wolff()
  use system
  implicit none
  integer :: i,j,k,s0,sc,cluster=0
  real(8) :: rn
  i=int(rn()*n)+1
  s0=spin(i)
  spin(i)=-spin(i)
  sc=i
  clusters=0
  do 
    do j=1,cor
      k=nbor(sc,j)
      if (k.ne.0) then
        if ((rn().lt.bp) .and. (spin(k).eq.s0)) then
            spin(k)=-spin(k)
            cluster=cluster+1
            clusters(cluster)=k
        endif
      endif
    enddo
    if(cluster==0) exit
    sc=clusters(cluster)
    cluster=cluster-1
  enddo
endsubroutine

subroutine cleardata
    use system
    arm=0;asm=0;m_4=0;binder=0;
    engy=0;engy2=0;
    test0val=0;
endsubroutine cleardata

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

