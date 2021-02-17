module otter_fibers_globals

    use otter_globals

    real(kind=DBL),allocatable,dimension(:,:) :: fibers
    real(kind=DBL),allocatable,dimension(:,:,:) :: contacts,fiber_path,control_points,bent_to_stright
    integer,allocatable,dimension(:) :: num_contacts,bent
    character(len=2) :: geocode='F1'

    contains



end module otter_fibers_globals