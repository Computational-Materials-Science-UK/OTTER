!!! OTTER - otter_fibers.f90
!
! OTTER: Complex Structure Generation Toolkit
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
! otter_spheres.f90 - module providing otter_fibers structure generation version
!	vb.1 - Working version
!   va.2 - Sean McDaniel
!   va.1 - Beck original with contributions from Skylar Mays

! External Dependencies:
!	None
! Internal Dependencies:
!	otter_main.f90

module fibers_place

    contains

    !!!!!
    ! otter_fibers - main subroutine for OTTER version dropping cylindrical fibers
    !	vb.1 - Working version
    !   va.2 - Sean McDaniel
    !   va.1 - Beck original with contributions from Skylar Mays

    ! External Dependencies:
    !	none
    ! Internal Dependencies:
    !	[TBA]
    subroutine otter_fibers(in_unit,out_unit)

        implicit none 

        !!!! Calling parameters
        integer,intent(IN) :: in_unit,out_unit  ! I/O unit numbers from otter_main

        !!!! Constants
        !!!! Variables

        !!!! Begin code...

        write(out_unit,*) ' In otter_fibers stub...'

    end subroutine otter_fibers

end module fibers_place

