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

        ! Requires otter globals, provided in module...
 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Variables & Constants


        if (debug) write(*,*) ' In otter_input:get_gen_input'

        write(out_unit,*) ' File management...'
        write(out_unit,'(a)',ADVANCE='NO') '   How many RVEs should be built for this set?   '
        read(in_unit,*) rves_num
        write(out_unit,'(a)',ADVANCE='NO') '   What is the batch name for this set of RVEs?  '
        read(in_unit,*) rves_batch_name
        write(out_unit,'(a)',ADVANCE='NO') '   Full path w/ trailing slash for this set?     '
        read(in_unit,'(a)') out_path
        write(out_unit,*) ' Universal settings...'
        write(out_unit,'(a)',ADVANCE='NO') '   Box dimensions (x, y, z):                     '
        read(in_unit,*) box_length(1:3)

        if (debug) write(*,*) ' Done otter_input:get_gen_input'

    end subroutine get_gen_input

    subroutine get_spheres_input(min_rad,max_rad,box_boundary)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Dependencies

        ! Requires otter globals, provided in module...
 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Variables & Constants
        real(kind=DBL)                  :: min_rad,max_rad
        real(kind=DBL), dimension(6)    :: box_boundary

        if (debug) write(out_unit,*) ' get_spheres_input(): enter'

        write(out_unit,*) ' Spheres specific inputs...'
        write(out_unit,'(a)',ADVANCE='NO') '   Minimum and maximum sphere radii?             '
        read(in_unit,*) min_rad, max_rad
        write(out_unit,'(a)',ADVANCE='NO') '   Boundary thicknesses (+/-x, +/-y, +/-z):      '
        read(in_unit,*) box_boundary(1:6)
        write(out_unit,'(a)',ADVANCE='NO') '   What is the min overlap distance for contact? '
        read(in_unit,*) min_olp
        write(out_unit,'(a)',ADVANCE='NO') '   What is the step size for moving primitives?  '
        read(in_unit,*) step_size

        if (debug) write(out_unit,*) ' get_spheres_input(): exit'

    end subroutine get_spheres_input

    subroutine get_fibers_input(min_rad, max_rad, min_length, max_length, friction)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Dependencies

        ! Requires otter globals, provided in module...
 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Variables & Constants
        real(kind=DBL)                  :: min_rad, max_rad, min_length, max_length, friction

        write(out_unit,*) ' Fiber specific inputs... '
        write(out_unit,'(a)',ADVANCE='NO') '   Minimum and maximum fiber radii?              '
        read(in_unit,*) min_rad, max_rad
        write(out_unit,'(a)',ADVANCE='NO') '   Minimum and maximum fiber length?             '
        read(in_unit,*) min_length, max_length
        write(out_unit,'(a)',ADVANCE='NO') '   What is the min overlap distance for contact? '
        read(in_unit,*) min_olp
        write(out_unit,'(a)',ADVANCE='NO') '   What is the step size for moving primitives?  '
        read(in_unit,*) step_size
        write(out_unit,'(a)',ADVANCE='NO') '   Min angle (deg) below horiz. for s/gliding?   '
        read(in_unit,*) friction

    end subroutine get_fibers_input

    subroutine get_ligaments_input()
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Dependencies

        ! Requires otter globals, provided in module...
        ! Requires otter_ligaments_globals: min_rad, max_rad, sphere_buffer, sphere_inflation, sphere_num, box_vol_factor, num_ligs
 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        use otter_ligaments_globals

        !! Variables & Constants

        !! Get Input
        write(out_unit,*) ' Ligaments specific inputs... '
        write(out_unit,'(a)',ADVANCE='NO') '   Number of spheres in compute box?   '
        read(in_unit,*) sphere_num
        write(out_unit,'(a)',ADVANCE='NO') '   (Nearest neighbor) connectivity number?   '
        read(in_unit,*) num_ligs
        write(out_unit,'(a)',ADVANCE='NO') '   ''Big Box'' volume factor (>1.0)?   '
        read(in_unit,*) box_vol_factor
        write(out_unit,'(a)',ADVANCE='NO') '   Sphere (node) min/max radii?   '
        read(in_unit,*) min_rad,max_rad
        write(out_unit,'(a)',ADVANCE='NO') '   Sphere (node) buffer and inflation?   '
        read(in_unit,*) sphere_buffer, sphere_inflation

    end subroutine get_ligaments_input


end module otter_input