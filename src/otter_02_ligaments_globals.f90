module otter_ligaments_globals

    use otter_globals

    real(kind=DBL),allocatable,dimension(:,:)     :: spheres,lig_map
    real(kind=DBL),allocatable,dimension(:,:,:)   :: nn,ligaments
    character(len=2) :: geocode='L1'

    real(kind=DBL)  :: min_rad,max_rad,sphere_buffer,sphere_inflation,box_vol_factor
    integer         :: num_ligs,sphere_num

    contains



end module otter_ligaments_globals