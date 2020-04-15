module otter_spheres_globals

    use otter_globals

    real(kind=DBL) :: min_rad, max_rad, sphere_buffer
    real(kind=DBL) :: shift, box_range
    real(kind=DBL),allocatable,dimension(:,:) :: spheres
    real(kind=DBL),allocatable,dimension(:,:,:) :: nn
    real(kind=DBL),allocatable,dimension(:,:) :: nn_count
    real(kind=DBL),allocatable,dimension(:,:) :: info

    integer :: sphere_num

    character :: lig_edge

    contains



end module otter_spheres_globals