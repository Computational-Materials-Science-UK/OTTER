!!! OTTER: Complex Structure Generation Toolkit
!	Maintained by the Computational Materials Science Research Group at the University of Kentucky, Dr. Matthew J. Beck, PI.
!	https://www.beckdt.engr.uky.edu
!
! OTTER master: Matthew Beck, m.beck@uky.edu
! OTTER developers: OTTER-dev team at GitHub.org
!
! https://github.com/Computational-Materials-Science-UK/OTTER
!
!!!

!!!!!
! otter_test.f90 - test program file for OTTER toolkit
!	vb.1 - Working version

! External Dependencies:
!	None
! Internal Dependencies:
!	otter_spheres.f90
!   otter_fibers.f90
!   otter_ligaments.f90

program otter_test

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Dependencies
    use otter_globals
    use spheres_place
    use fibers_place
    use ligaments_place

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Variables & Constants

    implicit none 

    ! constants
    

    ! variables
    integer :: option

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Begin code

    ! Initialize rand:
    write (*,*) ' In OTTER test, ',version
    call init_random()
    write(*,*) ' init_random() completed... ',rand()
    
    write(*,*) ''
    write(*,*) ' Please select a test option: '
    write(*,*) '    1. Test sphere rotation (w/ test contact)... '
    write(*,*) '    2. Test cone rotation (w/ test contact)... '
    read(*,*) option

    select case (option)
    case (1)
        call test_sphere_rotate()
    case(2)
        call test_cone_rotate()
    case default
        write(*,*) ' Invalid selection...'
    end select


end program otter_test

subroutine test_sphere_rotate()

    use otter_spheres_globals
    use spheres_place
    implicit none 

    integer :: sphere_num,sphere_max
    real(kind=DBL) :: step_size,min_olp

    sphere_max=2
    sphere_num=2
    step_size=0.1
    min_olp=0.01
    allocate(spheres(sphere_max,6))
    allocate(nn(sphere_num,5,3))

    spheres(:,:)=0.d0
    spheres(:,4)=1.d0
    !> test contact at 1,1,1 direction
    !spheres(2,1:3)=sqrt((1.98d0**2)/3.d0) 
    !> tests contact just above fall criteria in 0,1,0 direction
    spheres(2,3)=0.05d0
    spheres(2,2)=sqrt((1.98d0**2)-(0.05d0**2)) 

    debug=.true.
    call old_spheres(2,min_olp)
    write(*,*) ' Num contact after old_spheres: ',spheres(2,5)
    spheres(2,1:3)=rotate_sphere(2,step_size)
    write(*,*) ' New position after rotate_sphere: ',spheres(2,1:3), norm2(spheres(2,1:3))
    spheres(2,1:3)=rotate_sphere(2,step_size)
    write(*,*) ' New position after rotate_sphere: ',spheres(2,1:3), norm2(spheres(2,1:3))
    debug=.false.

end subroutine

subroutine test_cone_rotate()

    use otter_spheres_globals
    use spheres_place
    implicit none 

    integer :: sphere_num,sphere_max
    real(kind=DBL) :: step_size,min_olp

    sphere_max=3
    sphere_num=3
    step_size=0.1
    min_olp=0.01
    allocate(spheres(sphere_max,6))
    allocate(nn(sphere_num,5,3))

    spheres(:,:)=0.d0
    spheres(:,4)=1.d0
    !> Test sphere above and midway between two kissing spheres along Y
    spheres(2,2)=2.d0
    spheres(3,2)=1.d0
    spheres(3,3)=sqrt((1.98d0**2)-1.d0)
    !> contacted 1 @ origin, contacted 2 kisses 1 above in z and displace -0.25 in Y
    !> contacting is 1.98 from each of them, should roll off of the lower one, not cone rotate.
    spheres(2,2)=-0.25
    spheres(2,3)=sqrt(2.0**2-0.25**2)
    spheres(3,2)=1.56182
    spheres(3,3)=1.21701

    debug=.true.
    call old_spheres(3,min_olp)
    write(*,*) ' Num contact after old_spheres: ',spheres(3,5)
    spheres(3,1:3)=rotate_sphere(3,step_size)
    write(*,*) ' New position after rotate_sphere: ',spheres(3,1:3), norm2(spheres(3,1:3))
    call old_spheres(3,min_olp)
    write(*,*) ' Num contact after old_spheres AFTER 1st Rotation: ',spheres(3,5)
    spheres(3,1:3)=rotate_sphere(3,step_size)
    write(*,*) ' New position after rotate_sphere: ',spheres(3,1:3), norm2(spheres(3,1:3))
    debug=.false.

end subroutine
