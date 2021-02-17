!!! OTTER - otter_main.f90
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
! otter_main.f90 - main program file for OTTER toolkit
!	vb.1 - Working version

! External Dependencies:
!	None
! Internal Dependencies:
!	otter_spheres.f90
!   otter_fibers.f90
!   otter_ligaments.f90

program otter

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
    call init_random()

    ! Check for batch input, set-up input, output
    in_unit=5   ! Default for stdin
    out_unit=6  ! Default for stdout

    ! Batch input file must be named otter.in, batch output will be otter.out
    inquire(FILE='otter.in', exist=batch_exists)
    inquire(FILE='otter.out', exist=out_exists)

    if (batch_exists) then
        if (out_exists) then
            write(*,*) ' Batch input exists (otter.in). But otter.out exists and will be overwritten.'
            stop
        end if
        in_unit=51
        out_unit=52
        open(unit=in_unit,file='otter.in',status='OLD',action='READ')
        open(unit=out_unit,file='otter.out',status='NEW',action='WRITE')
        write(out_unit,*) ' Batch input/output selected.'
    end if

    if (debug) write(*,*) ' I/O setup: (in_unit, out_unit) ',in_unit,out_unit

    ! Write code header... 
    write(out_unit,*) ' Welcome to OTTER: Complex Structure Generation Toolkit v',VERSION
    write(out_unit,*) '   For questions, bugs, etc., visit: ',MASTER_WEB
    write(out_unit,*) '   OTTERmaster: ',MASTER_NAME,' ',MASTER_EMAIL
    write(out_unit,*) 

    ! Get OTTER structure type...
    write(out_unit,*) ' Available ''seed'' types: '
    write(out_unit,*) '     1. Dropped Spheres'
    write(out_unit,*) '     2. Dropped Fibers'
    write(out_unit,*) '     3. Ligamented/Swiss Cheese'
    !!!!! Add new options here, and below in select case...
    write(out_unit,*) ' Enter your selection: '
    read(in_unit,*) structure_type

    ! Call specific subroutine
    select case (structure_type)
    case (1)
        call otter_spheres()
    case (2)
        call otter_fibers()
    case (3)
        call otter_ligaments()
    !!!!! Add new options here to match above...
    case default
        write(out_unit,*) ' This option is invalid.'
        stop
    end select

    if (debug) write(*,*) ' Returned from otter_*'

end program
