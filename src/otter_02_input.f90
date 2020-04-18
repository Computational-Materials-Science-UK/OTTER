!!! OTTER - otter_input.f90
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
! otter_input.f90 - module providing general input handling to the OTTER codebase
!	vb.1 - Working version

! External Dependencies:
!	None
! Internal Dependencies:
!	otter_main.f90
!   otter_globals.mod

module otter_input

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Dependencies
    use otter_globals

    implicit none

    contains

    subroutine get_gen_input()

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Dependencies
 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Variables & Constants

        if (debug) write(*,*) ' In otter_input:get_gen_input'

        write (out_unit,*) ' How many RVEs do you want to build? '
        read(in_unit,*) rves_num
        write (out_unit,*) ' What is the batch name of your models? '
        read (in_unit,*) rves_batch_name
        write (out_unit,*) ' What is the full path with trailing slash for your models? '
        read(in_unit,'(a)') out_path
        write(out_unit,*) ' What is your desired box side length? '
        read(in_unit,*) box_length

        if (debug) write(*,*) ' Done otter_input:get_gen_input'

    end subroutine get_gen_input


end module otter_input