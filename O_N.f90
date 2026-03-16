module system
    implicit none

    integer,parameter :: o_n=5
    integer :: l                   
    integer :: n,step1,step2,step3
    integer, allocatable :: nbor(:,:)                         
    real(8), allocatable :: spin(:,:)
    real(8), allocatable :: m_dim(:),vect(:), clusters(:)
    real(8) :: temp,arm, asm, m_4,binder,c_r,c_r2,c_r3,c_r4
    real(8),parameter :: pi=acos(-1.d0)

end module system

program main
    use system
    implicit none
    integer :: i,j,k

    call initran(1)
    open(1,file='read.in',status='old')
    read(1,*)l,temp,step1,step2,step3
    close(1)

    call initial()
    do i=1,step1
        call wolff()
    enddo

    do i=1,step2
        call cleardata         
        do j=1,step3           
            call wolff()
            call sample   
        enddo                                                             
        call writers
    enddo

    call record()
    deallocate (spin)
    deallocate (nbor)
    deallocate (m_dim)
    deallocate (vect)
    deallocate (clusters)
end program

subroutine initial()
    use system
    implicit none
    integer i,j,x,y,z
    real(8) :: norm,rn
    logical :: spin_exist

    n=l**3
    allocate (nbor(0:n-1,6))
    allocate (spin(0:n-1,o_n))
    allocate (m_dim(o_n),vect(o_n))
    allocate (clusters(n))
    inquire(file='spin.dat',exist=spin_exist)
    if(spin_exist) then
        call readin()
    else
        do  i=0,n-1  
        98 continue
        do j=1,o_n
            spin(i,j)=2.d0*rn()-1.d0
        enddo
        norm=sum(spin(i,:)**2)
        if(norm>1.d0)goto 98
        spin(i,:)=spin(i,:)/sqrt(norm)
        enddo
    endif

    do i=0,n-1
        x=mod(i,l)
        y=mod(i/l,l)
        z=i/(l*l)
        nbor(i,1)=mod(x+1,l)+y*l+z*l*l
        nbor(i,2)=mod(x-1+l,l)+y*l+z*l*l
        nbor(i,3)=x+mod(y+1,l)*l+z*l*l
        nbor(i,4)=x+mod(y-1+l,l)*l+z*l*l
        nbor(i,5)=x+y*l+mod(z+1,l)*l*l
        nbor(i,6)=x+y*l+mod(z-1+l,l)*l*l
    end do
end subroutine initial

subroutine writers
    use system
    implicit none

    open(2,file='out.dat',position='append')    
    arm=arm/(dble(step3)*dble(n))
    asm=asm/(dble(step3)*dble(n)**2)
    m_4=m_4/(dble(step3)*dble(n)**4)
    binder=m_4/asm**2
    c_r=c_r*2.d0/(dble(step3)*dble(n))
    c_r2=c_r2*2.d0/(dble(step3)*dble(n))
    c_r3=c_r3*2.d0/(dble(step3)*dble(n))
    c_r4=c_r4*2.d0/(dble(step3)*dble(n))
    write(2,1)temp,arm,asm,m_4,binder,c_r,c_r2,c_r3,c_r4
    1 format(f5.3,8f18.12)
    close(2)
end subroutine writers

subroutine sample
    use system
    implicit none
    real(8) :: comp,m
    integer :: i,j
    integer :: x,y,z,x1,y1,z1
    m_dim=0;m=0;
    do i=0,n-1
        m_dim=m_dim+spin(i,:)
    enddo
    do i=0,n/2-1
        x=mod(i,l);y=mod(i/l,l);z=i/(l*l)
        x1=mod(x+l/2,l);y1=mod(y+l/2,l);z1=mod(z+l/2,l)
        comp=sum(spin(i,:)*spin(x1+y1*l+z1*l*l,:)) 
        c_r=c_r+comp
        if(comp> 1.d0) comp= 1.d0
        if(comp<-1.d0) comp=-1.d0
        c_r2=c_r2+cos(2.d0*acos(comp))
        c_r3=c_r3+cos(3.d0*acos(comp))
        c_r4=c_r4+cos(4.d0*acos(comp))
    enddo
    m=m+sum(m_dim**2)
    asm=asm+m 
    m_4=m_4+m**2 
    arm=arm+sqrt(m)  
end subroutine

subroutine cleardata
    use system
    arm=0;asm=0
    m_4=0
    binder=0
    c_r=0 
    c_r2=0            
    c_r3=0
    c_r4=0
endsubroutine cleardata

subroutine mcstep()
    use system
    implicit none
    integer :: i,j,k
    real(8) :: rn,ei,ef
    do i=0,n-1
        ef=0.d0
        ei=0.d0
    !       do
    !         call random_number(vect)
    !           vect=2.d0*vect-1.d0
    !           if (sum(vect**2) <= 1) exit
    !       enddo
        do k=1,o_n
            vect(k)=sqrt(-2.d0*log(rn()))*cos(2.d0*pi*rn()) !normal(0.d0,1.d0)
        enddo
        vect=vect/sqrt(sum(vect**2))
        do j=1,6
            ei=ei+sum(spin(i,:)*spin(nbor(i,j),:))
            ef=ef+sum(vect(:)*spin(nbor(i,j),:))
        enddo    
        if (rn()<=exp((ef-ei)/temp)) spin(i,:)=vect(:)
    enddo
endsubroutine

subroutine wolff()
    use system
    implicit none

    integer :: i,j,k,cluster=0,sc
    real(8) :: rn,k_c,k_n

    i=int(rn()*n)
    do k=1,o_n
        vect(k)=sqrt(-2.d0*log(rn()))*cos(2.d0*pi*rn()) !normal(0.d0,1.d0)
    enddo
    vect=vect/sqrt(sum(vect**2))
    k_c=sum(spin(i,:)*vect(:))
    spin(i,:)=spin(i,:)-2.d0*k_c*vect(:)
    clusters=0
    sc=i
    do 
        do j=1,6
        k_n=sum(spin(nbor(sc,j),:)*vect(:))
        k_c=sum(spin(sc,:)*vect(:))
        if (rn().lt.(1-exp(2.d0*k_c*k_n/temp))) then
            spin(nbor(sc,j),:)=spin(nbor(sc,j),:)-2.d0*k_n*vect(:)
            cluster=cluster+1
            clusters(cluster)=nbor(sc,j)
        endif
        enddo
        sc=clusters(cluster)
        if(cluster==0) exit
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

