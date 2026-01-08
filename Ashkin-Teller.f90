!----------------------!
 module systemvariables
   implicit none

   integer :: l,n,nn               
   integer,parameter :: z=4 !coordination number
   real(8) :: J2,temp
   real(8) :: J4
   integer, allocatable :: spin(:,:)     
   integer, allocatable :: istn(:)
   integer, allocatable :: neib(:,:)
   integer :: aaa(3)
   real(8) :: mmm(3),bp(-1:1)
   real(8) :: en1=0.d0
   real(8) :: en2=0.d0
   real(8) :: msigma=0.d0,msigma2=0.d0,msigma4=0.d0
   real(8) :: mtau=0.d0,mtau2=0.d0,mtau4=0.d0
   real(8) :: mst=0.d0,mst2=0.d0,mst4=0.d0
   integer :: bins,isteps,mcsteps,nmsr=0
   real(8) :: a1,a2,a3

 end module systemvariables
!--------------------------!

program AT
!----------------------------------!
   use systemvariables; implicit none

   integer :: i,j,k,nm,nw

   open(7,file='read.in',status='old')
   read(7,*)l,J2,J4,isteps,bins,mcsteps
   close(7)
   ! J2=1/temp
   nm=0;nw=1

   call initran(1)
   call initialize
   do i=1,isteps
      call mc(nm,nw)
   end do
   nmsr=0
   do j=1,bins
      do i=1,mcsteps
         call mc(nm,nw)
         call measure
      end do
      call writebindata()
      call checkpoint()
   end do
   call record()

   deallocate (neib,spin,istn)
end program AT
!-------------------!

 subroutine mc(nm,nw)
!----------------------------------!
 use systemvariables; implicit none
 integer :: nm,nw,st,ts
 integer :: nstack
 integer :: i,j,k,i1
 integer :: js,ks,inb,icsp
 real(8) :: dE
 real(8),external :: ran
!  do  i1=1,nm
   !   call metropolis
!  enddo
!  do i1=1,nw
     if(J2>J4)then
         call wolff1
     else
         call wolff2
     endif
!  enddo
 !-------------------------------
 contains
 !-------------------------------
 subroutine metropolis
 do k=1,n
   i=int(ran()*n)+1
   dE=0.d0
   st=int(2*ran())+1
   if (spin(st,i) == 0) cycle   
   do j=1,z
      dE=dE+2*(spin(st,i)*spin(st,neib(j,i)))*J2
      dE=dE+2*(spin(1,i)*spin(1,neib(j,i))*spin(2,i)*spin(2,neib(j,i)))*J4
   end do
   if(ran()<dexp(-dE)) spin(st,i)=-spin(st,i)
 enddo
 nmsr=nmsr+1
 end subroutine metropolis
 
 subroutine wolff1
   st=int(2*ran())+1
   ts=3-st
   i=int(ran()*n)+1
   do while (spin(st,i) == 0)
      i=int(ran()*n)+1
   end do
   nmsr=nmsr+1
   icsp=spin(st,i)
   spin(st,i)=-icsp
   nstack=1
   istn(1)=i
   do while (nstack > 0)
      js=istn(nstack)
      nstack=nstack-1
      do inb=1,z
         ks=neib(inb,js)
         if (spin(st,ks) == 0) cycle
         if((spin(st,ks) == icsp) .and. (ran() < bp(spin(ts,js)*spin(ts,ks)))) then
            spin(st,ks)=-icsp
            nstack=nstack+1
            istn(nstack)=ks
         end if
      end do
   end do
 end subroutine wolff1

 subroutine wolff2
   spin(3,:)=spin(1,:)*spin(2,:)
   st=int(2*ran())+1
   ts=3-st
   i=int(ran()*n)+1
   do while (spin(3,i) == 0)
      i=int(ran()*n)+1
   end do
   nmsr=nmsr+1
   icsp=spin(3,i)
   spin(3,i)=-icsp
   nstack=1
   istn(1)=i
   do while (nstack > 0)
      js=istn(nstack)
      nstack=nstack-1
      do inb=1,z
         ks=neib(inb,js)
         if (spin(3,ks) == 0) cycle
         if((spin(3,ks) == icsp) .and. (ran() < bp(spin(ts,js)*spin(ts,ks)))) then
            spin(3,ks)=-icsp
            nstack=nstack+1
            istn(nstack)=ks
         end if
      end do
   end do
   spin(st,:)=spin(3,:)*spin(ts,:)
 end subroutine wolff2
end subroutine mc
!---------------------!
    
 subroutine measure
!----------------------------------!
 use systemvariables; implicit none

 real(8) :: ms,mt,st
 integer :: i,j,k,e,s2,s4
 
 spin(3,:)=spin(1,:)*spin(2,:)
 s2=0;s4=0
 do i=1,n
   aaa=sum(spin(:,neib(:,i)),dim=2)
   s2=s2-spin(1,i)*aaa(1)-spin(2,i)*aaa(2)
   s4=s4-spin(3,i)*aaa(3)
 end do

 a1=(s2*J2+s4*J4)/(dble(nn)*2.d0)
 en1=en1+a1
 en2=en2+a1**2

 mmm=sum(spin(:,1:n),dim=2)/dble(nn)
 ms=mmm(1)
 msigma=msigma+dabs(ms)
 ms=ms**2
 msigma2=msigma2+ms
 msigma4=msigma4+ms**2

 mt=mmm(2)
 mtau=mtau+dabs(mt)
 mt=mt**2
 mtau2=mtau2+mt
 mtau4=mtau4+mt**2

 st=mmm(3)
 mst=mst+dabs(st)
 st=st**2
 mst2=mst2+st
 mst4=mst4+st**2

 end subroutine measure
!----------------------!

!---------------------------------!
subroutine writebindata()
 use systemvariables; implicit none

 en1=en1/dble(nmsr)
 en2=en2/dble(nmsr)

 msigma=msigma/dble(nmsr)
 msigma2=msigma2/dble(nmsr)
 msigma4=msigma4/dble(nmsr)
 msigma4=msigma4/msigma2**2
 msigma4=(1.d0-msigma4/3)

 mtau=mtau/dble(nmsr)
 mtau2=mtau2/dble(nmsr)
 mtau4=mtau4/dble(nmsr)
 mtau4=mtau4/mtau2**2
 mtau4=(1.d0-mtau4/3)*1.5d0

 mst=mst/dble(nmsr)
 mst2=mst2/dble(nmsr)
 mst4=mst4/dble(nmsr)
 mst4=mst4/mst2**2
 mst4=(1.d0-mst4/3)*1.5d0

 open(1,file='out.dat',position='append')
 write(1,1)l,nn,J2,J4,en1,en2,msigma,msigma2,msigma4,mtau,mtau2,mtau4,mst,mst2,mst4
 1 format(i6,'  ',i9,'  ',f8.5,'  ',f8.5,11f18.12)
 close(1)

 en1=0.d0;en2=0.d0
 msigma=0.d0
 msigma2=0.d0
 msigma4=0.d0
 mtau=0.d0
 mtau2=0.d0
 mtau4=0.d0
 mst=0.d0
 mst2=0.d0
 mst4=0.d0
 nmsr=0

end subroutine writebindata
!---------------------------!


subroutine initialize
!--------------------------!
   use systemvariables
   implicit none
   logical :: spin_exist
   integer :: i,j,k,ns,cnt
   integer, allocatable:: idx(:)
   real(8),external :: ran
   real(8) :: pp=0.75
   n=l**2
   nn=n
   allocate(idx(n))
   allocate(spin(1:3,0:n))
   allocate(istn(n))
   allocate(neib(1:z,1:n))
   inquire(file='spin.dat',exist=spin_exist)
   if(.false.) then
       call readin()
     else
         spin=0
         ns=nint(n*pp)
         do i=1,n
            idx(i)=i
         enddo
         do i=n,2,-1
            j=int(ran()*dble(i))+1
            k=idx(j)
            idx(j)=idx(i)
            idx(i)=k
         enddo
         do i=1,ns
            spin(1,idx(i))=2*int(2.d0*ran())-1
            spin(2,idx(i))=2*int(2.d0*ran())-1
         enddo
         nn=ns

         ! do i=1,n
         !    if (ran()>pp) then
         !       nn=nn-1
         !       cycle
         !    endif
         !    spin(1,i)=2*int(2.d0*ran())-1
         !    spin(2,i)=2*int(2.d0*ran())-1
         ! end do
   endif
   spin(1:3,0)=0
   bp(-1)=1-exp(-2.d0*abs(J2-J4))
   bp( 0)=0.d0
   bp( 1)=1-exp(-2.d0*(J2+J4))

   ! open(21,file='nbor.txt',status='old')
   ! do i=1,n
   !    read(21,*) neib(1,i),neib(2,i),neib(3,i),neib(4,i)
   !    if (all(neib(:,i).eq.0)) then
   !       spin(1:3,i)=0
   !       nn=nn-1
   !    end if
   ! end do
   ! close(21)
   
   cnt=1
   k=1
   do j=1,l
      do i=1,l
         neib(1,cnt)=(k-1)*l*l+(j-1)*l+mod(i,l)+1  !x+1,y,z
         neib(2,cnt)=(k-1)*l*l+(j-1)*l+l-mod(l-i+1,l) !x-1,y,z
         neib(3,cnt)=(k-1)*l*l+mod(j,l)*l+i    !x,y+1,z
         neib(4,cnt)=(k-1)*l*l+(l-mod(l-j+1,l)-1)*l+i  !x,y-1,z
         cnt=cnt+1
      end do
   end do 

   ! cnt=1
   ! do k=1,l
   !    do j=1,l
   !          do i=1,l
   !             neib(1,cnt)=(k-1)*l*l+(j-1)*l+mod(i,l)+1
   !             neib(2,cnt)=(k-1)*l*l+(j-1)*l+l-mod(l-i+1,l)
   !             neib(3,cnt)=(k-1)*l*l+mod(j,l)*l+i
   !             neib(4,cnt)=(k-1)*l*l+(l-mod(l-j+1,l)-1)*l+i
   !             neib(5,cnt)=mod(k,l)*l*l+(j-1)*l+i
   !             neib(6,cnt)=(l-mod(l-k+1,l)-1)*l*l+(j-1)*l+i
   !             cnt=cnt+1
   !          end do
   !    end do
   ! end do

end subroutine initialize
!-------------------------!

subroutine record()
  use systemvariables
  open(12,file='spin.dat')
  write(12,*)spin
  close(12)
end subroutine record
subroutine readin()
  use systemvariables
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

!----------------------!
 real(8) function ran()
!----------------------------------------------!
! 64-bit congruental generator                 !
! iran64=oran64*2862933555777941757+1013904243 !
!----------------------------------------------!
 implicit none

 real(8)    :: dmu64
 integer(8) :: ran64,mul64,add64
 common/bran64/dmu64,ran64,mul64,add64

 ran64=ran64*mul64+add64
 ran=0.5d0+dmu64*dble(ran64)

 end function ran

!---------------------!
 subroutine initran(w)
!---------------------!
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
 
!----------------------!
