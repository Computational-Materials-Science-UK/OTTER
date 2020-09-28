module otter_globals

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Global constants...
    integer,parameter :: DBL=SELECTED_REAL_KIND(7)
    character(len=3),parameter :: VERSION='b.1'
    character(len=*),parameter :: MASTER_NAME='Matthew Beck'
    character(len=*),parameter :: MASTER_EMAIL='m.beck@uky.edu'
    character(len=*),parameter :: MASTER_WEB='https://beckdt.engr.uky.edu'
    real(kind=DBL),parameter :: PI=4.d0*atan(1.d0)
    
    logical :: DEBUG=.false.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Global variables...
    integer                         :: in_unit,out_unit
    integer                         :: rves_num
    character(len=24)               :: rves_batch_name
    character(len=192)              :: out_path
    real(kind=DBL)                  :: min_olp, step_size
    real(kind=DBL),dimension(3)     :: box_length
    logical,parameter               :: test_random=.TRUE.       !!! SET THIS TO TRUE to get same rand numbers!

    contains

    !!! OTTER - subroutine init_random
    !   Generates an 8 digit integer seed using date_and_time, then seeds a call to srand.
    !   Executed once at the beginning of OTTER.

    !	vb.1 - Working version

    ! Requires:
    !   None

    subroutine init_random()

        integer,dimension(8)            :: values       ! stores results from date_and_time
        integer                         :: seed_int     ! 8 digit int for seed srand...

        call date_and_time(values=values)
        seed_int=values(8)*100000+values(7)*1000+values(6)*10+values(5)
        if (debug) write(*,*) ' init_random: values(1:8), seed_int: ',values,seed_int

        if (test_random) then
            write(*,*) ' *****WARNING*****'
            write(*,*) ' You are running the test_random version.  See otter_globals!'
            write(*,'(a)',ADVANCE='NO') ' Enter 8-digit random seed: '
            read(*,*) seed_int
        end if
        call srand(seed_int)

        if (debug) write(*,*) ' init_random: seed_int: ',seed_int

    end subroutine init_random

end module