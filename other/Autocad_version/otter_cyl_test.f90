program otter_cyl_test

    use otter_cyl_rotate
    implicit none

    integer :: test_num=0
    real(kind=8),dimension(3)           :: ellipse_center
    real(kind=8)                       :: B,A  ! minor/major axes of ellipse
    real(kind=8),dimension(3)           :: r    ! point through which plane passes

    write(*,*) ' What do you want to test?'
    write(*,*) '   1. get_R'
    write(*,*) '   2. get_ellipse'
    write(*,*) '   3. get_t (calls get_ellipse3!)'
    write(*,*) '   4. rotate_beck (calls all)'
    write(*,*) '   5. other tests'
    write(*,*) '   6. get_slope'
    write(*,*) ' Enter integer, any other number will exit: '
    read(*,*) test_num

    select case (test_num)
        case (1)
            CALL test_get_R !! WORKS
        case (2)
!            CALL test_get_ellipse(r,ellipse_center,B,A)  !! WORKS
        case (3)
            CALL test_get_t(r,ellipse_center,B,A)
        case(4)
            CALL test_get_rotate
        case(5)
            CALL test_other
        case(6)
            CALL test_get_slope
    end select

    write (*,*) ' Exiting...'

end program

subroutine test_other
    real(kind=8),dimension(3) :: v=(/1.d0,1.d0,1.d0/)

    write(*,*) ' vec: ',v(1:3)
    write(*,*) ' norm2: ',norm2(v(1:3))

end subroutine

subroutine test_get_slope
    use otter_cyl_rotate
    implicit none

    real(kind=8)                :: t,phi,a,b,slope,curvature
    real(kind=8),dimension(3)   :: ellipse_center
    real(kind=8),parameter      :: PI=4.0*atan(1.d0)

    t=-0.1
    phi=0.d0
    a=1.d0
    b=1.5d0
    ellipse_center=(/0.d0,0.d0,0.d0/)
    CALL get_ellipse_slope(t,phi,ellipse_center,a,b,slope,curvature)
    write(*,*) ' 1 - slope,curvature: ',slope,curvature
    t=1.001*PI
    CALL get_ellipse_slope(t,phi,ellipse_center,a,b,slope,curvature)
    write(*,*) ' 2 - slope,curvature: ',slope,curvature

end subroutine

subroutine test_get_rotate
    use otter_cyl_rotate
    implicit none

    real(kind=8),dimension(3,3)      :: q
    real(kind=8),dimension(3,3)   :: p
    real(kind=8),dimension(2,3)     :: rp
    real(kind=8),dimension(3)        :: r,qq,p_uvec
    real(kind=8)                     :: pmag,qmag


    p(1,:)=(/0.,0.,1.9/)
    p(2,:)=(/0.,1.,2./)
    p(3,:)=p(2,:)-p(1,:)
    pmag=norm2(p(3,:))
    p_uvec=p(3,:)/pmag
    r=(/0.501d0,0.499d0,1.d0/)
    rp(1,:)=p(1,:)-r(:)
    rp(2,:)=p(2,:)-r(:)
    q(1,:)=(/1.,1.,0./)
    q(2,:)=(/0.,0.,0./)
    q(3,:)=q(2,:)-q(1,:)
    qmag=norm2(q(3,:))
    qq=dot_product((r-q(1,:)),q(3,:))/dot_product(q(3,:),q(3,:))*q(3,:)+q(1,:)

    CALL rotate_beck(r,qq,p,q,rp,p_uvec,pmag,qmag)

end subroutine

subroutine test_get_t(r,ellipse_center,B,A)
    use otter_cyl_rotate
    implicit none

    !calling parameters
    real(kind=8),dimension(3),intent(inout)           :: ellipse_center
    real(kind=8),intent(inout)                        :: B,A  ! minor/major axes of ellipse
    real(kind=8),dimension(3),intent(inout)           :: r    ! point through which plane passes

    !local parameters
    real(kind=8)                            :: t,p_0,z_0,phi

    CALL test_get_ellipse(r,ellipse_center,B,A,phi,t)
    write(*,*) 'ellipse_center, r: ',ellipse_center,r
    write(*,*) 'B,A: ',B,A
    !CALL get_t(r,ellipse_center,B,A,t,p_0,z_0)

    write(*,*) 't,phi: ',t,phi

end subroutine

subroutine test_get_ellipse(r,ellipse_center,B,A,phi,t)
    use otter_cyl_rotate
    implicit none

    !calling parameters
    real(kind=8),dimension(3),intent(inout)           :: ellipse_center
    real(kind=8),intent(inout)                        :: B,A  ! minor/major axes of ellipse
    real(kind=8),dimension(3),intent(inout)           :: r    ! point through which plane passes

    !local parameters
    real(kind=8),dimension(3)           :: vec  ! with z-axis defines plane orientation
    real(kind=8),dimension(3,3)         :: q    ! cylinder: q(1,:)-start pt, q(2,:)-end pt, q(3,:)-vector centerline
    real(kind=8)                        :: qmag,slope,curvature ! length of cylinder
    real(kind=8),dimension(3)           :: qq   ! point on centerline from which a radius hits r
    real(kind=8),intent(inout)                        :: phi,t


    vec=(/1.,0.,0./)
    r=(/0.2d0,9.d0,-1.d0/)
    q(1,:)=(/1.,1.,0./)
    q(2,:)=(/0.,0.,0./)
    q(3,:)=q(2,:)-q(1,:)
    qmag=norm2(q(3,:))
    qq=dot_product((r-q(1,:)),q(3,:))/dot_product(q(3,:),q(3,:))*q(3,:)+q(1,:)
    write(*,*) 'qq: ',qq,precision(qq(1))

    CALL get_ellipse3(vec,r,q,qmag,qq,ellipse_center,A,B,phi,t,slope,curvature)

    CALL get_ellipse_slope(t,phi,ellipse_center,A,B,slope,curvature)

!    write(*,*) 'ellipse_center: ',ellipse_center
!    write(*,*) 'B,A: ',B,A

end subroutine

subroutine test_get_R
    use otter_cyl_rotate
    implicit none

    real(kind=8),dimension(3)       :: a,b,c
    real(kind=8),dimension(3,3)     :: RRR
    a=(/-0.41771793314805167,       0.20972647051203022,        7.3537081238809641E-004/)
    b=(/-0.89367306125317691,       0.44869248374328358,        4.8491878349628060E-003/)
    a=(/1.d0, 0.d0, 0d0/)
    b=(/0.d0,1.d0,0.d0/)
    c=(/2.d0,0.d0,1.d0/)
    CALL get_R(a,b,RRR)

    write(*,*) 'a,b,matmul(RRR,a): ',a,b,matmul(RRR,a)
    write(*,*) 'c,matmul(RRR,c): ',c,matmul(RRR,c)

end subroutine
