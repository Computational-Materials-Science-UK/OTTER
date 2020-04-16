module otter_globals

    integer,parameter :: DBL=SELECTED_REAL_KIND(7)
    character(len=3),parameter :: VERSION='b.1'
    character(len=*),parameter :: MASTER_NAME='Matthew Beck'
    character(len=*),parameter :: MASTER_EMAIL='m.beck@uky.edu'
    character(len=*),parameter :: MASTER_WEB='https://beckdt.engr.uky.edu'
    
    logical,parameter :: DEBUG=.false.

    contains

end module