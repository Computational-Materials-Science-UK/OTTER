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
    logical :: batch_exists,out_exists  ! checks for existing batch I/O files
    integer :: structure_type           ! User selected structure version

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Begin code

    ! Initialize rand:
    write (*,*) ' In OTTER test, ',version

    call init_random()
    write (*,*) ' init_random() completed... ',rand()


end program otter_test