!!! OTTER - otter_ligaments.f90
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
! otter_spheres.f90 - module providing otter_ligaments structure generation version
!	vb.1 - Working version
!   va.1 - Beck with contributions from Mujan Seif
!   v0.1 - Original OTTER variant, with thanks to many...

! External Dependencies:
!	None
! Internal Dependencies:
!	otter_main.f90

module ligaments_place

    contains

    !!!!!
    ! otter_ligaments - main subroutine for OTTER version dropping cylindrical ligaments
    !	vb.1 - Working version
    !   va.1 - Beck with contributions from Mujan Seif
    !   v0.1 - Original OTTER variant, with thanks to many...

    ! External Dependencies:
    !	none
    ! Internal Dependencies:
    !	[TBA]
    subroutine otter_ligaments(in_unit,out_unit)

        implicit none 

        !!!! Calling parameters
        integer,intent(IN) :: in_unit,out_unit  ! I/O unit numbers from otter_main

        !!!! Constants
        !!!! Variables

        !!!! Begin code...

        write(out_unit,*) ' In otter_ligaments stub...'
        write(*,*) 'Learning git'

    end subroutine otter_ligaments

end module ligaments_place