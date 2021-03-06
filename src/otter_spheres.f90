!!! OTTER - otter_spheres.f90
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
! otter_spheres.f90 - module providing otter_spheres structure generation version
!	vb.1 - Working version
!   va.1 - Based on Skylar Mays work, with additions by Katie Moody

! External Dependencies:
!	None
! Internal Dependencies:
!	otter_main.f90
!   otter_input.mod

module spheres_place

    use otter_globals
    use otter_input
    use otter_spheres_globals

    contains

    !!! OTTER - subroutine otter_spheres
    !   Spheres with a range of sizes (min_rad, max_rad) fall into a 'big box' that is box_boundary
    !   larger than a cube of box_length. Spheres contact each other when they overlap by between
    !   min_olp and max_olp, and rotate around each contacted sphere until they contact 3 spheres,
    !   at which point they are fixed. Spheres fall when they rotate off of another sphere.

    !	vb.1 - Working version
    !   va.1 - Based on Skylar Mays work, with additions by Katie Moody

    ! Requires:
    !   otter_input
    !   otter_spheres_globals
    !   otter_globals

    subroutine otter_spheres()

        implicit none

        ! local variables
        real(kind=DBL)                  ::  min_rad, max_rad, &
                                            big_box_vol,sphere_range, &
                                            nn_avg,sphere_num_in
        real(kind=DBL),dimension(6)     ::  box_boundary,bigbox
        real(kind=DBL),dimension(3)     ::  box_range, shift
        integer                         ::  i,j, &
                                            sphere_num, excess_num, &
                                            scr_unit,dat_unit,rve_scr_unit,rve_nnc_unit, &
                                            rve,num_bad

        character(len=200)              ::  full_scr_name,full_nnc_name,full_name, &
                                            full_rves_batch_name,full_nncd_name, &
                                            saveas,exportstl
        character(len=50)               ::  date_time
        character(len=8)                ::  num_char
        logical                         ::  in_box,in_bounds




        ! version history
        ! 1: random spheres with buffer, falland connect without rotating, nn count
        ! 2: attempt to rotate until touches 3 sheres
        ! 3: implentation of subroutines
        ! b.1: Beck re-write

        write (out_unit,*) ' Welcome to Spheres vb.1!'
        write (out_unit,*) ' '

        ! Get Global/General input, see otter_input
        call get_gen_input()
        ! Get spheres input..., see otter_input
        call get_spheres_input(min_rad,max_rad,box_boundary)

        ! compute 'big box' lengths, volume; shift origin of cut box to (0,0,0)
        big_box_vol=1.d0
        do i=1,3
            box_range(i)=box_length(i)+box_boundary(i*2-1)+box_boundary(i*2)
            big_box_vol=big_box_vol*box_range(i)
            shift(i)=0.d0-box_boundary(i*2)
        end do
        ! Estimate max number of spheres as big_box_vol/smallest sphere
        sphere_num=int(big_box_vol/( 4.d0*PI*(min_rad-min_olp)**3/3.d0 ))
        excess_num=int( (( sqrt( (box_length(3)/2/max_rad)-3 ) +3 )*2*max_rad)* &
                   box_length(1)*box_length(2) / ( 4.d0*PI*(min_rad-min_olp)**3/3.d0 ) )
        write(out_unit,*) ''
        write(out_unit,*) ' Big box lengths, volume: ',box_range(:),big_box_vol
        if (debug) write(out_unit,*) ' (debug) shift, sphere_num: ',shift(:),sphere_num

        ! Generate structures...
        write (out_unit,*) ' Generating Structure(s)...'

        !!! Allocate arrays:
        ! spheres:  sphere data for each of sphere_num spheres - nn_count is number of contacting spheres...
        !       x-pos, y-pos, z-pos, radius, nn_count, in_box (1=yes, 0=no)
        allocate(spheres(sphere_num,6))

        ! nn:  nearest neighbor table for each of sphere_num shperes:
        !   contact sphere number, distance, overlap distance
        ! Each sphere in theory has 3 nn bc it falls until it touches 3..., 
        ! NN are only detailed here for the sphere landing, though all contacts are counted in spheres(i,5)
        allocate(nn(sphere_num,5,3))

        !info:  volume and nn data for all RVEs...
        !	avgnn, total volume overlap
        allocate(info(rves_num,2))
        info(:,:)=0.d0

        ! Set-up output for Auto-CAD multi-script
        ! **** THERE is a weakness here... output should be centralized, and it should check for existing files before crashing upon discovery....
        scr_unit=53 ! was 52
        dat_unit=54 ! was 47?
        full_scr_name=trim(rves_batch_name)//".scr"
        full_nnc_name=trim(rves_batch_name)//".out"
        open(unit=scr_unit,file=trim(full_scr_name),status='new',action='write')
        open(unit=dat_unit,file=trim(full_nnc_name),status='new',action='write')

        ! Write dat file header...
        call fdate(date_time)
        write(dat_unit,*) rves_batch_name
        write(dat_unit,*) ' Date generated: ',date_time
        write(dat_unit,*) ' Number of RVEs: ',rves_num
        write(dat_unit,*) ' Box length, Boundaries (+x, -x, +y, -y, +z, -z):'
        write(dat_unit,*) box_length(:),box_boundary(:)
        write(dat_unit,*) ' Min rad, max rad, min overlap, step_size:'
        write(dat_unit,*) min_rad, max_rad, min_olp, step_size
        write(dat_unit,*) ''
        write(dat_unit,*) ' RVE data...'

        !Begin building structures, one model at a time...e
        do rve=1,rves_num
            write(out_unit,*) 'RVE #',rve,'...'
            
            ! Set-up output for individual RVE
            write(num_char,"(i4)") rve !cast a four character version of the rve number
            full_name=trim(rves_batch_name)//trim(adjustl(num_char))
            full_rves_batch_name=trim(rves_batch_name)//trim(adjustl(num_char))//".out"
            full_nncd_name=trim(rves_batch_name)//trim(adjustl(num_char))//".plain"
            rve_scr_unit=63 ! was 50
            rve_nnc_unit=64 ! was 46
            open(unit=rve_scr_unit,file=trim(full_rves_batch_name),status='new',action='write')
            open(unit=rve_nnc_unit,file=trim(full_nncd_name),status='new',action='write')

            !initialize...
            spheres(:,:)=0.d0
            spheres(:,6)=1.d0
            nn(:,:,:)=0.d0
            i=1
            num_bad=0

            sphere_range=max_rad-min_rad
            if(debug) write(*,*) 'otter_spheres: sphere_range: ',sphere_range

            ! Main loop... create spheres until sphere maximum is exceeded or until the last 0.1*sphere_num have been bad, that is, above the big_box.
            do while ((i .le. sphere_num) .and. (num_bad .le. excess_num))
                if (debug) write(out_unit,'(i3,a4)',ADVANCE='NO') i,'... '

                ! randomly place a sphere above box
                do j=1,2
                    spheres(i,j)=rand()*box_range(j)+shift(j)
                end do
                spheres(i,3)=box_range(3)+2*max_rad+shift(3)
                spheres(i,4)=min_rad+rand()*sphere_range
                if (debug) write(out_unit,*) spheres(i,:)

                in_bounds=.true.

                ! Drop sphere until contact, or hits bottom.
                do while ((spheres(i,3) .gt. 0.d0+shift(3)) .and. (spheres(i,5) .lt. 3) .and. &
                  in_bounds)

                    ! Check for contact
                    CALL old_spheres (i)
                    if (debug) write(*,*) 'Re-entered program spheres_place from old_spheres'

                    select case (int(spheres(i,5)))
                    case (0)
                        ! Drop sphere by step_size
                        spheres(i,3)=spheres(i,3)-step_size
                    case (1)
                        spheres(i,1:3)=rotate_sphere(i)
                    case (2)
                        spheres(i,1:3)=rotate_sphere(i)
                    case default
                        if (debug) write(out_unit,*) ' Stopped...'
                    end select

                    ! Check if off an edge...
                    do j=1,2
                        if ( (spheres(i,j) .lt. shift(j)) .or. (spheres(i,j) .gt. &
                          box_range(j)+shift(j)) ) in_bounds=.false.
                    end do

                end do ! sphere done falling and rotating...

                !write(*,'(a,i5,6f10.4)') ' Sphere done, i, pos, nn, in_box:',i,spheres(i,:)
                !read(*,*)

                ! Check if sphere is in final box...
                do j=1,3
                    if ((spheres(i,j) .lt. (0.d0-spheres(i,4))) .or. (spheres(i,j) .gt. &
                      box_length(j)+spheres(i,4))) then
                        spheres(i,6) = 0.d0
                        if (debug) write(*,*) ' Out of box: i, pos: ',i,spheres(i,1:3)
                    end if
                end do

                ! Check for sphere still above box...
                if (spheres(i,3) .gt. box_range(3)+shift(3)) then
                    ! increment count of consecutive above box spheres, re-initialize sphere i for re-use...
                    num_bad=num_bad+1
                    spheres(i,:)=0.d0
                    spheres(i,6)=1.d0
                    nn(i,:,:)=0.d0
                else
                    num_bad=0
                    i=i+1
                end if

            end do ! Done creating spheres for RVE
            if (debug) write(out_unit,*) 'generating spheres complete'

            ! Check if sphere_num was too low...
            if (num_bad .le. excess_num) then
                write(out_unit,*) ' RAN OUT of SPHERES!  Sphere_num too low, and probably other problems.'
                stop
            end if

            ! Compute NN statistics
            sphere_num_in=0
            nn_avg=0
            do j=1,i-1
                if (spheres(j,6) .gt. 0.5d0) then
                    sphere_num_in=sphere_num_in+1
                    nn_avg=nn_avg+spheres(j,5)
                end if
            end do
            nn_avg=nn_avg/sphere_num_in

            !!!!! Output RVE information....
            ! NOTE: Overlap data is IGNORED!
            ! NN data is also not written!!!

            write(dat_unit,*) ' RVE #',rve,': ',full_name, ' - Avg NN: ',nn_avg

            ! Write RVE to batch scr file...
            ! Write and union spheres...
            202 format('_sphere ',f0.3,',',f0.3,',',f0.3,' ',f0.3)
            write(scr_unit,'(a)') '-osnap off'
            write(scr_unit,'(a)') '_zoom -300,-300,-300 500,500,500'
            do j=1,i-1
                if (spheres(j,6) .gt. 0.5d0) write(scr_unit,202) spheres(j,1:4)
            end do
            write(scr_unit,'(a)') '_union all '
            write(scr_unit,'(a)') '_group create network  all '

            ! Calculate and write bigbox...
            do j=1,3
                bigbox(j)=0.d0+shift(j)-2.d0*max_rad
                bigbox(j+3)=box_range(j)+shift(j)+2.d0*max_rad
            end do
            301 format('_box ',f0.4,',',f0.4,',',f0.4,' ',f0.4,',',f0.4,',',f0.4)
            write(scr_unit,301) bigbox(1:6)
            write(scr_unit,'(a)') '_group create bigbox  last '

            ! Calculation box...
            write(scr_unit,301) 0.,0.,0.,box_length(1),box_length(2),box_length(3)
            write(scr_unit,'(a)') '_group create smallbox  last '

            ! Subtract groups, generate cut structure, move to first quadrant
            write(scr_unit,'(a)') '_subtract g bigbox  g smallbox '
            write(scr_unit,'(a)') '_subtract g network  g bigbox '
            write(scr_unit,'(a)') '_move all  0,0,0 1,1,1'

            ! Save data
            saveas='_saveas  '//trim(out_path)//trim(adjustl(full_name))//'.dwg'
            exportstl='_export '//trim(out_path)//trim(adjustl(full_name))//'.stl all'
            write(scr_unit,'(a)') 'FILEDIA 0'
            write(scr_unit,'(a)') trim(saveas)
            write(scr_unit,'(a,a)') trim(exportstl),' '

            ! Clean up...
            write(scr_unit,'(a)') '_erase all '
            write(scr_unit,'(a)') '-purge g network y y'
            write(scr_unit,'(a)') '-purge g bigbox y y'
            write(scr_unit,'(a)') '-purge g smallbox y y'

        end do  ! models

        ! Clean up RVE scr...
        write(scr_unit,'(a)') 'FILEDIA 1'
        close(unit=scr_unit)

    end subroutine otter_spheres

    !!! OTTER - function rotate_sphere - returns new_pos
    !   Finds axis of rotation:  For 1 contacting sphere this is the normal to the plane containing both z_hat and the vector from the contacted to the contacting sphere.  For 2 contacting spheres, this is vector between the two contacted spheres.
    !   Rotates point at center of contacting sphere about axis by a distance that depends the distance of the center of the contacting sphere to the axis.

    !	vb.1 - Working version
    !   va.1 - Based on Skylar Mays work, with additions by Katie Moody

    ! Requires:
    !   otter_spheres_globals
    !   otter_globals

    function rotate_sphere (i) result(new_pos)
        use otter_spheres_globals
        use otter_globals
        use otter_math
        implicit none

        ! Declare calling parameters:
        integer, intent(in) :: i

        !Declare local variables
        real(kind=DBL),dimension(3)     :: new_pos,pt,axis,npt,temp_pos
        real(kind=DBL)                  :: theta,dist
        integer                         :: j

        if(debug) write(out_unit,*) 'Sphere SUBROUTINE: ',i,spheres(i,:)
        
        select case (int(spheres(i,5)))
        case (1)
            if (debug) write (*,*) ' rotate_sphere: in case 1'

            new_pos=do_sphere_rotation(spheres(i,1:3),spheres(int(nn(i,1,1)),1:3),nn(i,1,2))

        case (2)
            if (debug) write (*,*) ' rotate_sphere: in case 2'

            ! pt & axis set origin @ 1st contacted sphere
            ! axis is vector between 1st/2nd contacted sphere
            pt=spheres(i,1:3)-spheres(int(nn(i,1,1)),1:3)
            axis=spheres(int(nn(i,2,1)),1:3)-spheres(int(nn(i,1,1)),1:3)
            ! Find dist to axis for theta rotation calc...
            call nearpt_pt2line3(spheres(i,1:3),spheres(int(nn(i,2,1)),1:3),spheres(int(nn(i,1,1)),1:3),npt,dist)
            theta=step_size/dist

            if (debug) then
                write(*,*) ' rotate_sphere 2: axis: ',axis
                write(*,*) ' rotate_sphere 2: theta: ',theta
                write(*,*) ' rotate_sphere 2: pt: ',pt,norm2(pt)
            end if

            ! axis in ax2om must be unit vector to avoid expansion/contraction during rotation!
            axis=axis/norm2(axis)
            new_pos=matmul(ax2om_d([axis(1),axis(2),axis(3),theta]),pt)
            new_pos=new_pos+spheres(int(nn(i,1,1)),1:3)

            do j=1,2
                temp_pos=do_sphere_rotation(spheres(i,1:3),spheres(int(nn(i,j,1)),1:3),nn(i,j,2))
                if ( ( norm2(temp_pos-spheres( int( nn(i,mod(j,2)+1,1) ),1:3 )) .ge. nn(i, &
                  mod(j,2)+1,2) ) .and. ( temp_pos(3) .lt. new_pos(3) ) ) then
                    if (debug) write(*,*) ' Better option is: j: ',j
                    new_pos=temp_pos
                end if
            end do
            
        case default
            if (debug) write (*,*) ' rotate_sphere: in case default'
            write(out_unit,*) ' In rotate_sphere: Too many contacts (spheres(i,5))=',spheres(i,5)
        end select

        if (debug) write(*,*) ' Distance moved: ',norm2(new_pos-spheres(i,1:3))
        

        ! Test to make sure rotation lowers sphere...
        !! There is a small error in this test... it is possible that for contacts near the
        !! top of the sphere, this could go the wrong way...
        if (new_pos(3) .gt. spheres(i,3)) then
            theta=0.d0-theta
            new_pos=matmul(ax2om_d([axis(1),axis(2),axis(3),theta]),pt)+spheres(int(nn(i,1,1)),1:3)
        end if

    END function rotate_sphere

    FUNCTION do_sphere_rotation(pos_contacting,pos_contacted,dist) result(new_pos)
        use otter_spheres_globals
        use otter_globals
        use otter_math
        implicit none 

        ! Declare calling parameters
        real(kind=DBL),intent(in),dimension(3)    :: pos_contacting,pos_contacted
        real(kind=DBL),intent(in)                 :: dist

        ! Variables
        real(kind=DBL),dimension(3)               :: new_pos,axis,pt
        real(kind=DBL)                            :: theta

        ! pt sets origin at contacted sphere
        pt=pos_contacting-pos_contacted

        ! check if next to or below...
        if (pt(3) .le. 0.d0) then
            if (debug) write(*,*) ' do_sphere_rotation: sphere falling...'
            new_pos=pt-[0.d0,0.d0,step_size]+pos_contacted
        else ! axis is pt vector cross z-hat, set theta to step / separation distance...
            axis=cross3(pt,[0.d0,0.d0,1.d0])
            theta=step_size/dist

            if (debug) then
                write(*,*) ' do_sphere_rotate: axis: ',axis
                write(*,*) ' do_sphere_rotate: theta: ',theta
                write(*,*) ' do_sphere_rotate: pt: ',pt,norm2(pt)
            end if

            ! Check for sphere on top...
            if (norm2(axis) .eq. 0.d0) axis=cross3([1.d0*rand(),1.d0*rand(),pt(3)],[0.d0,0.d0,1.d0])
            ! Comput rotated position...
            ! axis in ax2om must be unit vector to avoid expansion/contraction during rotation!
            axis=axis/norm2(axis)
            new_pos=matmul(ax2om_d([axis(1),axis(2),axis(3),theta]),pt)+pos_contacted
        end if

    END FUNCTION do_sphere_rotation


    SUBROUTINE nn_detailed (dist,i,j)
        !Purpose to save information of the detailed nn list
        use otter_spheres_globals
        use otter_globals
        implicit none

        !declare calling parameters
        integer,intent(in) :: i,j
        real(kind=DBL),dimension(4),intent(in) :: dist

        if (spheres(i,5) .gt. 5) then
            write(out_unit,*) ' TOO MANY SPHERE NN''S for NN table!'
        else
            nn(i,int(spheres(i,5)),1) = j
            nn(i,int(spheres(i,5)),2) = dist(4)
            nn(i,int(spheres(i,5)),3) = overlap_volume(dist,i,j)
        end if

        if (debug) write (*,*) 'Detailed nnc list used: ', nn(i,:,:)

        RETURN
    END SUBROUTINE nn_detailed

    SUBROUTINE old_spheres (i)
        !..Purpose: called in to check new sphere against old sphere
        use otter_spheres_globals
        use otter_globals
        implicit none

        !Declare calling parameters
        integer, intent(in) :: i

        !Declare local variables
        real(kind=DBL), dimension(4) :: dist
        integer :: j

        if (debug) write(*,*) 'Old Sphere SUBROUTINE'

        ! Reset spheres to find new set of contacts!
        spheres(i,5)=0
        j=1

        ! Check if sphere is +x/-x/+y/-y out of

        do while (j .lt. i)
            !calc distances between sphere (i) and sphere (j)
            !The vector dist points from old sphere j to new sphere i
            dist(1:3)=spheres(i,1:3)-spheres(j,1:3)
            dist(4)=norm2(dist(1:3))
            if (debug) write(out_unit,*) 'Old_spheres: NN loop: ',j,dist(4)

            !Check and see if sphere separation is more than min_olp
            if (dist(4)-spheres(i,4)-spheres(j,4) .le. (0.d0-min_olp)) then
                spheres(i,5)=spheres(i,5)+1
                spheres(j,5)=spheres(j,5)+1
                if (debug) write (out_unit,*) 'Old_spheres: New NN (j,i,spheres(j,5),spheres(i,5)): ',j,i,spheres(j,5),spheres(i,5)

                !save information to detailed nn chart
                CALL nn_detailed (dist,i,j)

            end if

            j=j+1
        end do ! j .lt. i

        if (debug) then
            write(*,*) 'End of Old Sphere SUBROUTINE'
            write(*,*) ' nn:',nn(i,int(spheres(i,5)),:)
        end if

        RETURN

    END SUBROUTINE old_spheres

    real(kind=DBL) FUNCTION overlap_volume (dist,i,j)
        !...Purpose: called to caluclate overlap volume between two spheres
        !!!!*****THIS function may be unnecessary, and is wrong because the h fraction for each sphere is wrong
        use otter_spheres_globals
        use otter_globals
    
        implicit none

        !Declare calling parameters
        integer, intent(in) :: i,j
        real(kind=DBL),dimension(4),intent(in) :: dist

        !Declare local variables
        real(kind=DBL),dimension(2) :: volume, h

        h(1)=(dist(4)-spheres(i,4)-spheres(j,4))*spheres(j,4)/(spheres(i,4)+spheres(j,4))
        h(2)=(dist(4)-spheres(i,4)-spheres(j,4))*spheres(i,4)/(spheres(i,4)+spheres(j,4))

        !Dome math of both spheres
        volume(1) = (1.d0/3.d0)*pi*(h(1)**2)*((3*spheres(i,4))-h(1))
        volume(2) = (1.d0/3.d0)*pi*(h(2)**2)*((3*spheres(j,4))-h(2))
        if (debug) write(*,*) 'Overlap: volume from i, j: ',volume
        overlap_volume = volume(1)+volume(2)

        RETURN
    end function overlap_volume


end module spheres_place