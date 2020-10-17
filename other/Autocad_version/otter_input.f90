module otter_input

implicit none

contains
    subroutine get_input(geocode,num_models,out_name,path,box,min_rad,max_rad,buffer,min_length,max_length,cylinder_num,lig_edge)

        !passed parameters
        character(len=10),intent(in) :: geocode
        character(len=200),intent(out) :: path
        character(len=18),intent(out) :: out_name
        character(len=1),intent(out) :: lig_edge
        integer,intent(out) :: num_models,cylinder_num
        real(kind=8),intent(out) :: box,min_rad,max_rad,buffer,min_length,max_length

        !local parameters
        !none...

        write(*,*) ''
        write(*,*) ' How many models do you want to make? '
        read(*,*) num_models
        write(*,*) ' What is the batch name of your models? '
        read(*,*) out_name
        write(*,*) ' What is the full path with trailing slash for your models? '
        read(*,'(a)') path
        write(*,*) ' What is your desired box side length? '
        read(*,*) box
        ! cylinder specific inquiry
        write(*,*) ' Please enter minimum and maximum cylinder radius: '
        read(*,*) min_rad, max_rad
        write(*,*) ' Please enter the max overlap distance: '
        read(*,*) buffer

        write(*,*) ' Please enter minimum and maximum cylinder length: '
        read(*,*) min_length, max_length
        write(*,*) ' Please enter the number of cylinders: '
        read(*,*) cylinder_num

        write(*,*) ' Ligamented edges? '
        read (*,*) lig_edge
    end subroutine

end module
