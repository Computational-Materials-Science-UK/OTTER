!!! OTTER - otter_fibers.f90
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
! otter_fibers.f90 - module providing otter_fibers structure generation version
!	vb.1 - Working version
!   va.2 - Sean McDaniel
!   va.1 - Beck original with contributions from Skylar Mays

! External Dependencies:
!	None
! Internal Dependencies:
!	otter_main.f90

module fibers_place

    use otter_globals
    use otter_input
    use otter_fibers_globals
    use otter_math

    contains

    !!!!!
    ! otter_fibers - main subroutine for OTTER version dropping cylindrical fibers
    !	vb.1 - Working version
    !   va.2 - Sean McDaniel
    !   va.1 - Beck original with contributions from Skylar Mays

    ! External Dependencies:
    !	none
    ! Internal Dependencies:
    !	[TBA]
    subroutine otter_fibers()

        implicit none 

        !!!! Calling parameters

        !!!! Constants
        !!!! Variables

        real(kind=DBL)                  ::  min_rad, max_rad, &
                                            big_box_vol,radius_range,length_range, &
                                            dir_rand, min_length, max_length, height, &
                                            max_height, friction
        real(kind=DBL),dimension(6)     ::  box_boundary,bigbox
        real(kind=DBL),dimension(3)     ::  box_range, shift
        integer                         ::  i,j,k, periodic_added, periodic_j, periodic_i, &
                                            max_fibers, go_glide, &
                                            scr_unit,dat_unit,rve_scr_unit,rve_nnc_unit, &
                                            rve,num_bad,excess_num,glide_num

        character(len=200)              ::  full_scr_name,full_nnc_name,full_name, &
                                            full_rves_batch_name,full_nncd_name, &
                                            saveas,exportstl
        character(len=50)               ::  date_time
        character(len=8)                ::  num_char
        logical                         ::  stopped, in_box, glide_in
        logical,dimension(2)            ::  stop_glide
        logical,parameter               ::  test_glide=.TRUE. ! Set TRUE to test... must be false for real structures!

        !!!! Begin code...

        write (out_unit,*) ' Welcome to Fibers vb.1!'
        write (out_unit,*) ' '

        if (test_glide) then
            write(*,*) ' *****WARNING*****'
            write(*,*) ' You are running the test_glide version, see otter_fibers'
            read(*,*)
        end if

        ! Get Global/General input, see otter_input
        call get_gen_input()

        ! Get fiber specific input...
        call get_fibers_input(min_rad, max_rad, min_length, max_length, friction, glide_num)

        ! Adjust inputs...
        glide_in=.false.
        if (glide_num .eq. 1) glide_in=.true.
        friction=asin(friction/360.d0*2.d0*PI)
        if (debug) write (out_unit,*) ' Friction in absolute z: ',friction
        box_range(1:3)=box_length(1:3)

        ! compute 'big box' lengths, volume; shift origin of cut box to (0,0,0)
        big_box_vol=1.d0
        ! don't need shift any more, awaiting removal...
        shift=0.d0
        do i=1,3
            big_box_vol=big_box_vol*box_range(i)
        end do
        ! Estimate max number of fibers as 4*big_box_vol/cylindrical plate containing smallest cylinder... no idea if this is any good...
        ! Adding 25% to account for basal boundary amount.
        max_fibers=5*int(big_box_vol/( PI*(min_length/2.)**2.*(2.*min_rad) ))
        ! set up for periodic images
        max_fibers=9*max_fibers
        write(out_unit,*) ''
        write(out_unit,*) ' Big box lengths, volume: ',box_range(:),big_box_vol
        if (debug) write(out_unit,*) ' (debug) shift, max_fibers: ',shift(:),max_fibers

        ! Generate structures...
        write (out_unit,*) ' Generating Structure(s)...'

        !!! Allocate arrays:
        ! spheres:  sphere data for each of max_fibers spheres - nn_count is number of contacting spheres...
        !       end1x, end1y, end1z, rad1, end2x, end2y, end2z, radz, length
        allocate(fibers(max_fibers,9))

        ! nn:  nearest neighbor table for each of max_fibers fibers:
        !   contact fiber number, contact_pt3D (x,y,z), separation distance
        ! Each fiber has 2 (3 with rolling) plus those on top contacts..., 
        allocate(contacts(max_fibers,7,5))
        allocate(num_contacts(max_fibers))

        radius_range=max_rad-min_rad
        length_range=max_length-min_length
        if(debug) write(*,*) 'ranges: ',radius_range,length_range

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
        write(dat_unit,*) box_length,box_boundary(:)
        write(dat_unit,*) ' Min rad, max rad, min_length, max_length, min overlap, step_size:'
        write(dat_unit,*) min_rad, max_rad, min_length, max_length, min_olp, step_size
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
            fibers(:,:)=0.d0
            contacts(:,:,:)=0.d0
            num_contacts(:)=0
            i=1
            num_bad=0
            max_height=box_range(3)

            radius_range=max_rad-min_rad
            length_range=max_length-min_length
            if(debug) write(*,*) 'otter_fibers: radius_range,length_range: ',radius_range,length_range

            ! Main loop... create fibers until fiber maximum is exceeded or until the last 0.05*max_fibers have been bad, that is, above the big_box.
            excess_num=int(0.05*max_fibers/9.d0)
            !!!!!!! TEMP ADJUST
            excess_num=1000000
            if(debug) write(out_unit,*) 'excess_num: ', excess_num
            do while ((i .lt. max_fibers) .and. (num_bad .le. excess_num))
                if (.true.) write(out_unit,'(i5,a4)',ADVANCE='NO') i,'... '

                ! randomly place a fiber above box
                fibers(i,1)=rand()*box_range(1)+shift(1)
                fibers(i,2)=rand()*box_range(2)+shift(2)
                fibers(i,3)=max_height+2.*max_rad+shift(3)
                fibers(i,4)=min_rad+rand()*radius_range
                fibers(i,9)=min_length+rand()*length_range
                dir_rand=rand()*PI
                fibers(i,5)=fibers(i,1)+fibers(i,9)*cos(dir_rand)
                fibers(i,6)=fibers(i,2)+fibers(i,9)*sin(dir_rand)
                fibers(i,7)=fibers(i,3)
                fibers(i,8)=min_rad+rand()*radius_range

                ! Drop fiber until contact, or hits bottom.
                stopped=.false.
                ! Enabling gliding for this fiber.
                stop_glide(:)=(.not. glide_in)
                if (debug) write(out_unit,*) 'stop_glide: ',stop_glide(:)
                do while (.not. stopped)

                    ! Check for contacts
                    call fiber_contact(i,max_rad)
                    if (debug) write(*,*) 'Re-entered program spheres_place from fiber_contact'

                    select case (num_contacts(i))
                    case (0)
                        ! Drop fiber by step_size
                        fibers(i,3)=fibers(i,3)-step_size
                        fibers(i,7)=fibers(i,7)-step_size
                    case (1)
                        fibers(i,:)=fiber_rotate(i,step_size)
                    case (2)
                        if (.not. (stop_glide(1) .and. stop_glide(2))) then
                            if (.true.) write(out_unit,*) ' Gliding!...'
                            fibers(i,:)=fiber_glide(i,step_size,friction,stop_glide)
                            if (debug) write(out_unit,*) ' Returned from gliding w/ stop_glide: ',stop_glide

                            ! IF test glide is enabled, stop after first glided fiber...
                            if (stop_glide(1) .and. stop_glide(2)) then
                                stopped=.true.
                                if (.true.) write(out_unit,*) ' Stopped with two contacts, both below friction. '
                            else
                                if (.true.) write(out_unit,*) ' Gliding again...'
                                if (test_glide) then
                                    write(out_unit,*) 'test_glide max_fiber invoked!',i
                                    max_fibers=i+8
                                    glide_in=.false.
                                end if
                            end if
                        else
                            if (.true.) write(out_unit,*) ' Stopped with two contacts and gliding turned off!'
                            stopped=.true.
                        end if

                    case default
                        ! Kludge here... in theory, there could be three contacts and fiber could still roll... this should be rare...
                        if (.true.) write(out_unit,*) ' Stopped w/ 3 or more contacts...',contacts(i,1:3,5)
                        stopped=.true.
                    end select

                    ! Check if reached bottom
                    !bottom=0.d0-0.25*box_length(3)
                    !if ((fibers(i,3) .le. bottom) .or. (fibers(i,7) .le. bottom)) then
                    !    stopped=.true.
                    !    write(*,*) 'Stopped by touching bottom...'
                    !end if

                end do ! dropping loop

                ! check if fiber is above box...
                height=max(fibers(i,3),fibers(i,7))
                if (height .gt. max_height) max_height=height
                if (min(fibers(i,3),fibers(i,7)) .gt. box_range(3)) num_bad=num_bad+1

                ! check fiber length...
                if (norm2(fibers(i,1:3)-fibers(i,5:7))-fibers(i,9) .gt. 0.001*max_length) then
                    write(*,*) ' FIBER LENGTH error! ',i
                    read(*,*)
                end if
                
                ! set up periodic images...
                if (debug) write(*,*) ' Entering periodic...'
                periodic_added=0
                do periodic_i=1,3
                    do periodic_j=1,3
                        if (.not. ((periodic_i .eq. 2) .and. (periodic_j .eq. 2))) then
                            i=i+1
                            periodic_added=periodic_added+1
                            fibers(i,:)=fibers(i-periodic_added,:)
                            fibers(i,1)=fibers(i,1)+(periodic_i-2)*box_range(1)
                            fibers(i,5)=fibers(i,5)+(periodic_i-2)*box_range(1)
                            fibers(i,2)=fibers(i,2)+(periodic_j-2)*box_range(2)
                            fibers(i,6)=fibers(i,6)+(periodic_j-2)*box_range(2)
                        end if
                    end do
                end do

                i=i+1

            end do ! fiber generation loop

            if (debug) write(out_unit,*) 'generating spheres complete'

            ! Check if sphere_num was too low...
            if (num_bad .le. excess_num) then
                write(out_unit,*) ' WARNING: excess_num not well selected - not a critical error.'
                !stop
            end if

            ! Write RVE to batch scr file...
            ! Write and union spheres...
            
            202 format('_sphere ',f0.3,',',f0.3,',',f0.3,' ',f0.3)
            203 format('_cone ',f0.3,',',f0.3,',',f0.3,' ',f0.3,' T ',f0.3,' A ',f0.3,',',f0.3,',',f0.3)
            write(scr_unit,'(a)') '-osnap off'
            write(scr_unit,'(a)') '_zoom -300,-300,-300 500,500,500'

            do j=1,i-1

                ! Check whether fiber has a least one end in box...
                in_box=.true.
                do k=1,3
                    if ( (min(fibers(j,k),fibers(j,k+4)) .gt. box_range(k)) .or. &
                         (max(fibers(j,k),fibers(j,k+4)) .lt. 0.d0) ) in_box=.false.
                end do

                ! If in the box, write to script...
                if (in_box) write(scr_unit,203) fibers(j,1:4),fibers(j,8),fibers(j,5:7)
            end do

            write(scr_unit,'(a)') '_union all '
            write(scr_unit,'(a)') '_group create network  all '

            ! Calculate and write bigbox...
            do j=1,3
                bigbox(j)=0.d0+shift(j)-2.d0*max_length
                bigbox(j+3)=box_range(j)+shift(j)+2.d0*max_length
            end do
            301 format('_box ',f0.4,',',f0.4,',',f0.4,' ',f0.4,',',f0.4,',',f0.4)
            write(scr_unit,301) bigbox(1:6)
            write(scr_unit,'(a)') '_group create bigbox  last '

            ! Calculation box...
            write(scr_unit,301) 0.,0.,0.,box_range(1:3)
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

            if (test_glide) call init_random()

        end do ! rves loop

        ! Clean up RVE scr...
        write(scr_unit,'(a)') 'FILEDIA 1'
        close(unit=scr_unit)
                

    end subroutine otter_fibers

    function fiber_rotate(i,step_size) result(new_pos)
        implicit none

        ! Declare calling parameters:
        integer, intent(in)             :: i
        real(kind=DBL),intent(in)       :: step_size

        !Declare local variables
        real(kind=DBL),dimension(3)     :: axis,fiber_dir
        real(kind=DBL),dimension(9)     :: new_pos
        real(kind=DBL),dimension(2,3)   :: end_pt,new_end_pt
        real(kind=DBL),dimension(2)     :: length
        real(kind=DBL),dimension(3,3)   :: R,Rold
        real(kind=DBL)                  :: theta
        integer                         :: direction
        logical                         :: drop_out

        if(debug) write(out_unit,*) 'fiber_rotate SUBROUTINE: ',i,fibers(i,:)
        ! Check for vertical...
        if (norm2(fibers(i,1:2)-fibers(i,5:6)) .le. 4.d0*step_size) then
            !DROP
            new_end_pt(1,1:2)=fibers(i,1:2)
            new_end_pt(2,1:2)=fibers(i,5:6)
            new_end_pt(1,3)=fibers(i,3)-step_size
            new_end_pt(2,3)=fibers(i,7)-step_size

            ! If this is turn on (drop_out=.true.), then vertical fibers are just dropped out of box.... warning is printed to the screen.
            drop_out=.true.
            if (drop_out) then
                new_end_pt(1,3)=new_end_pt(1,3) - 100000
                new_end_pt(2,3)=new_end_pt(2,3) - 100000
                write(*,*) ' Dropping OUT!'
            end if
        else

            ! We will rotate the vectors from the contact pt to the fiber axis end points
            end_pt(1,1:3)=fibers(i,1:3)-contacts(i,1,2:4)
            end_pt(2,1:3)=fibers(i,5:7)-contacts(i,1,2:4)

            ! Kludge, use distance from contact on surface to center line end as length
            length(1)=norm2(end_pt(1,1:3))
            length(2)=norm2(end_pt(2,1:3))

            ! Get the theta and axis of rotation...
            ! The axis is the vec perp to the plane containing the contacting fiber axis & (0,0,1)
            axis(:)=cross3([end_pt(2,1:3)-end_pt(1,1:3)],[0.d0,0.d0,1.d0])

            ! Set rotation direction and max step towards longer end...
            direction=2
            if (length(1) .ge. length(2)) direction=1
            theta=step_size/length(direction)

            ! Check if contact point is near the middle, then shift... Warning printed...
            fiber_dir(1:3)=fibers(i,1+2*direction:3+2*direction)/fibers(i,9)
            if (abs(length(1)-length(2)) .lt. 2.d0*step_size) then
                end_pt(1,1:3)=end_pt(1,1:3)+2.d0*step_size*fiber_dir(1:3)
                end_pt(2,1:3)=end_pt(2,1:3)+2.d0*step_size*fiber_dir(1:3)
                if (.true.) write (*,*) ' Shifting due to middle contact: ',fiber_dir
            end if

            !write(*,*) 'Length (long, short): ', direction, length(direction), length(mod(direction,2)+1)
            !write(*,*) 'Fibers end1:         ',fibers(i,1:3)
            !write(*,*) 'Fibers end2:         ',fibers(i,5:7)              
            !write(*,*) 'Contact i, contacts: ',i, contacts(i,1,1:5)
            !write(*,*) 'end_pt(1,1:3):       ',end_pt(1,1:3)
            !write(*,*) 'end_pt(2,1:3):       ',end_pt(2,1:3)

            ! normalize axis so ax2om keeps vector length constant...
            axis(:)=axis(:)/norm2(axis(:))
            !write(*,*) 'axis, theta:                ',axis,theta,norm2(axis(:))

            ! Compute rotation matrix
            R=ax2om_d([axis(1:3),theta])

            ! move end pt that should go down, check, reverse if necessary, move other pt
            new_end_pt(direction,:)=matmul(R,end_pt(direction,:))+contacts(i,1,2:4)
            if (new_end_pt(direction,3) .ge. end_pt(direction,3)+contacts(i,1,4)) then
                ! Reverse rotation by inverted rotation axis to keep theta positive
                axis=[0.d0,0.d0,0.d0]-axis
                if (debug) write(*,*) 'new theta: ',theta
                Rold=R
                R=ax2om_d([axis(1:3),theta])
                new_end_pt(direction,:)=matmul(R,end_pt(direction,:))+contacts(i,1,2:4)
            end if
            if (new_end_pt(direction,3) .ge. end_pt(direction,3)+contacts(i,1,4)) then
                write(*,*) ' ROTATE ERROR! No rotation lowers long end.'
                write(*,*) '   Likely fatal, quit code and notify OTTERMaster.',theta
                !write(*,*) R 
                !write(*,*) Rold
                write(*,*) '   Fiber new, old: ',new_end_pt(direction,:),end_pt(direction,:)+contacts(i,1,2:4)
                read(*,*)
            end if
            new_end_pt(mod(direction,2)+1,:)=matmul(R,end_pt(mod(direction,2)+1,:))+contacts(i,1,2:4)
            !write(*,*) ' Fiber new, old: ',new_end_pt(direction,:),end_pt(direction,:)+contacts(i,1,2:4)

            ! This is more "sliding" over contact then just rotating, address contact point drift...
            ! This is not the actual drift just an estimate... fiber is sliding towards longer end.
            new_end_pt(1,:)=new_end_pt(1,:)+step_size*fiber_dir(1:3)
            new_end_pt(2,:)=new_end_pt(2,:)+step_size*fiber_dir(1:3)
            
            ! Check for broken-ness in finding/incorporating contact point.
            if (length(1)+length(2) .gt. 55) then
                write(*,*) ' ROTATE ERROR! Length doesn''match up...'
                write(*,*) '   Likely fatal, quit code and notify OTTERMaster.'
                write(*,*) ' end_pts, contact: ',end_pt(1,1:3)
                write(*,*) '                   ',end_pt(2,1:3)
                write(*,*) '                   ',contacts(i,1,2:4)
                read(*,*)
            end if

        end if

        new_pos(1:3)=new_end_pt(1,1:3)
        new_pos(4)=fibers(i,4)
        new_pos(5:7)=new_end_pt(2,1:3)
        new_pos(8:9)=fibers(i,8:9)

    end function fiber_rotate


    subroutine fiber_contact(i,max_rad)
        implicit none 

        !Declare calling parameters
        integer, intent(in) :: i
        real(kind=DBL),intent(in) :: max_rad

        !Declare local variables
        real(kind=DBL), dimension(3)    :: contact_vec,contact_pt3D
        real(kind=DBL), dimension(2)    :: radius,near_ptP,sign_contact, near_endP
        real(kind=DBL), dimension(2,3)  :: near_pt3D,fiber_vec,face_vec, &
                                           center, near_end3D
        real(kind=DBL)                  :: dist,bottom, dist_end, sign_fiber
        integer                         :: j,k,contact_case,l
        logical                         :: contact,alt_contact=.false.

        num_contacts(i)=0
        contacts(i,:,:)=0

        if(debug) write(*,*) 'In fiber_contact...'

        ! Check gross contact
        do j=1,i-1
            contact=.true.

            do k=1,3
                if ( ( max(fibers(i,k),fibers(i,k+4)) .lt. &
                       min(fibers(j,k),fibers(j,k+4))-2.d0*max_rad ) .or. &
                     ( min(fibers(i,k),fibers(i,k+4)) .gt. &
                       max(fibers(j,k),fibers(j,k+4))+2.*max_rad ) ) contact=.false.
            end do

            ! check fine contact...
            if (contact) then
                if (debug) write (*,*) ' Check Fine contact: ',i,j
                contact=.false.
                
                call segments_dist_3d(fibers(i,1:3),fibers(i,5:7),fibers(j,1:3),fibers(j,5:7), &
                                      dist,near_pt3D(1,:), near_pt3D(2,:), near_ptP(1), near_ptP(2))

                radius(1)=fibers(i,4)+near_ptP(1)*(fibers(i,8)-fibers(i,4))
                radius(2)=fibers(j,4)+near_ptP(2)*(fibers(j,8)-fibers(j,4))
                center(1,:)=fibers(i,1:3)+near_ptP(1)*(fibers(i,5:7)-fibers(i,1:3))
                center(2,:)=fibers(j,1:3)+near_ptP(2)*(fibers(j,5:7)-fibers(j,1:3))
                contact_vec(:)=near_pt3D(2,:)-near_pt3D(1,:)
                fiber_vec(1,:)=(fibers(i,5:7)-fibers(i,1:3))/fibers(i,9)
                fiber_vec(2,:)=(fibers(j,5:7)-fibers(j,1:3))/fibers(j,9)

                contact_case=0
                if ((near_ptP(1) .eq. 0.d0) .or. (near_ptP(1) .eq. 1.d0)) contact_case=1
                if ((near_ptP(2) .eq. 0.d0) .or. (near_ptP(2) .eq. 1.d0)) contact_case=contact_case+2
                contact_case=contact_case+10

                sign_contact(1)=1.d0
                sign_contact(2)=-1.d0

                ! rounded end approximation...
                alt_contact=.false.
                if (.not. alt_contact) then
                    select case (contact_case)
                    case (3)
                        !! END-TO-END
                        ! Project contact_vec into each face to for face_vec with length
                        ! radius minus 1/2 min_olp -- this is a kludge with min_olp...
                        ! If these face_vec segments intersect (min dist=0) then contact.
                        ! Math for projection is vec-plane_vec*dot(vec,plane_vec)/norm(plane_vec)
                        
                        do k=1,2
                            sign_fiber=1.d0
                            if (near_ptP(k) .eq. 0.d0) sign_fiber=-1.d0

                            face_vec(k,1:3)=( sign_contact(k)*contact_vec(:) ) - &
                                            ( sign_fiber*fiber_vec(k,:) ) * &
                                            dot_product(sign_contact(k)*contact_vec(:), &
                                                        sign_fiber*fiber_vec(k,:)) / &
                                            norm2(sign_fiber*fiber_vec(k,:))
                            face_vec(k,:)=face_vec(k,:)/norm2(face_vec(k,:))*radius(k)-0.5*min_olp
                        end do
                        call segments_dist_3d(center(1,:),face_vec(1,:),center(2,:),face_vec(2,:), &
                                              dist_end, near_end3D(1,:), near_end3D(2,:), &
                                              near_endP(1), near_endP(2) )
                        if (dist_end .lt. 1.e-6) then
                            contact_pt3D(1:3)=( near_end3D(1,:)+near_end3D(2,:) ) / 2.d0
                            contact=.true. 
                        end if

                    case (2,1)
                        k=contact_case

                        sign_fiber=1.d0
                        if (near_ptP(k) .eq. 0.d0) sign_fiber=-1.d0

                        face_vec(k,1:3)=( sign_contact(k)*contact_vec(:) ) - &
                                        ( sign_fiber*fiber_vec(k,:) ) * &
                                        dot_product(sign_contact(k)*contact_vec(:), &
                                                    sign_fiber*fiber_vec(k,:)) / &
                                        norm2(sign_fiber*fiber_vec(k,:))
                        face_vec(k,:)=face_vec(k,:)/norm2(face_vec(k,:))*radius(k)
                        
                        l=i
                        if (k .eq. 2) l=j
                        call segments_dist_3d(center(k,:),face_vec(k,:),fibers(l,1:3), &
                                              fibers(l,5:7),dist_end,near_end3D(1,:), &
                                              near_end3D(2,:), near_endP(1), near_endP(2) )
                        if (dist_end .lt. radius(mod(k,2)+1)-min_olp) then
                            contact_pt3D(1:3)=near_end3D(2,:)
                            contact=.true. 
                        end if

                    case default

                        !!!! ORIGINAL ROUNDED END
                        if (dist .le. (radius(1)+radius(2)-min_olp)) then 
                            contact=.true.
                            contact_pt3D(1:3)=near_pt3D(1,1:3)+0.5*(near_pt3D(2,1:3)-near_pt3D(1,1:3))
                        end if
                        !!!! END ORIGINAL ROUNDED END

                    end select


                    
                else

                !near_pt3D(1,:)=fibers(i,1:3)+near_ptP(1)*(fibers(i,5:7)-fibers(i,1:3))
                !near_pt3D(2,:)=fibers(j,1:3)+near_ptP(2)*(fibers(j,5:7)-fibers(j,1:3))
                contact_vec(:)=near_pt3D(2,:)-near_pt3D(1,:)
                fiber_vec(1,:)=near_pt3D(1,:)-(fibers(i,5:7)+near_ptP(1)*(fibers(i,1:3)-fibers(i,5:7)))
                fiber_vec(2,:)=near_pt3D(2,:)-(fibers(j,5:7)+near_ptP(2)*(fibers(j,1:3)-fibers(j,5:7)))

                if ((near_ptP(1) .eq. 0.d0) .or. (near_ptP(1) .eq. 1.d0)) then
                    face_vec(1,:)=contact_vec(:)-dotn(contact_vec(:),fiber_vec(1,:)/norm2(fiber_vec(1,:)))
                    if ((near_ptP(2) .eq. 0.d0) .or. (near_ptP(2) .eq. 1.d0)) then
                        ! end-to-end : Check if proj of contact in each face is less than radii...
                        if (debug) write (*,*) ' checking end-to-end'
                        face_vec(2,:) = ([0.d0,0.d0,0.d0] - contact_vec(:)) - &
                                         dotn(([0.d0,0.d0,0.d0] - contact_vec(:)), &
                                         fiber_vec(2,:)/norm2(fiber_vec(2,:)))
                        ! kludge here... overlap is more than min_olp becuase it's applied to both face_vecs...
                        if ((norm2(face_vec(1,:)) .lt. radius(1)-min_olp) .and. (norm2(face_vec(2,:)) .lt. radius(2)-min_olp)) then
                            ! end-to-end contact
                            if (debug) write (*,*) ' end-to-end'
                            contact=.true.
                            contact_pt3D(:)=0.5*contact_vec(:)-dotn(0.5*contact_vec(:), &
                                            fiber_vec(1,:)/norm2(fiber_vec(1,:)))+near_pt3D(1,:)
                        end if
                    else
                        ! end-to-side : Check if edge point of i is within radius of j.
                        ! There is a small kludge here... the contact point is set to the corner of the contacting cylinder...
                        if (debug) write (*,*) ' checking end-to-side'
                        contact_pt3D(:)=near_pt3D(1,:)+radius(1)*face_vec(1,:)/norm2(face_vec(1,:))
                        if (pt2line(fibers(j,1:3),fibers(j,5:7),contact_pt3D(:)) .lt. radius(2)-min_olp) then
                            ! end-to-side contact
                            if (debug) write (*,*) ' end-to-side'
                            contact=.true.
                        endif
                    end if

                else 
                    if ((near_ptP(2) .eq. 0.d0) .or. (near_ptP(2) .eq. 1.d0)) then
                        ! end-to-side : Check if edge point of j is within radius of i.
                        ! There is a small kludge here... the contact point is set to the corner of the contacted cylinder...
                        contact_pt3D(:)=near_pt3D(2,:)+radius(2)*face_vec(2,:)/norm2(face_vec(2,:))
                        if (pt2line(fibers(i,1:3),fibers(i,5:7),contact_pt3D(:)) .lt. radius(1)-min_olp) then
                            ! end-to-side contact
                            contact=.true.
                        end if
                    else
                        ! side-to-side contact
                        if (dist .lt. radius(1)+radius(2)-min_olp) then
                            ! side-to-side contact
                            contact=.true.
                            contact_pt3D(:)=near_pt3D(1,:)+contact_vec(:)/dist*(radius(1)-0.5*(radius(1)+radius(2)-dist))
                        end if
                    end if
                end if
                end if

            end if ! end fine check

            if (contact) then
                if (.true.) then
                    write(*,*) ' Record fine contact!, Contact Case ',i,j,contact_case-10
                end if
                !write(*,*) ' Record fine contact! ',i,j
                !write(*,*) ' Fiber 1: ',fibers(i,1:3),fibers(i,5:7)
                !write(*,*) ' Fiber 2: ',fibers(j,1:3),fibers(j,5:7)
                !write(*,*) ' dist, near_ptP(1),(2): ',dist, near_ptP(1), near_ptP(2)
                !write(*,*) ' near_pt3D(1),(2): ',near_pt3D(1,:), near_pt3D(2,:) 
                num_contacts(i)=num_contacts(i)+1
                contacts(i,num_contacts(i),1)=j
                contacts(i,num_contacts(i),2:4)=contact_pt3D(1:3)
                contacts(i,num_contacts(i),5)=dist
            end if

        end do

        bottom=0.d0-0.25*box_length(3)
        !!!!!!!!  SHOULD ONLY BE FOR TEST_GLIDE....
        bottom=0
        do k=3,7,4
            if (fibers(i,k) .le. bottom) then
                num_contacts(i)=num_contacts(i)+1
                contacts(i,num_contacts(i),1)=0
                contacts(i,num_contacts(i),2:4)=fibers(i,k-2:k)
                contacts(i,num_contacts(i),5)=0.d0
                write(*,*) ' One end is touching the bottom...', k
            end if
        end do

    end subroutine

    function fiber_glide(i,step_size, friction, stop_glide) result(new_pos)
        implicit none

        ! Declare calling parameters:
        integer, intent(in)             :: i
        real(kind=DBL),intent(in)       :: step_size,friction

        ! Output
        real(kind=DBL),dimension(9)     :: new_pos
        logical,dimension(2)            :: stop_glide

        ! Local variables...
        integer                         :: j,one
        real(kind=DBL),dimension(2,3)   :: top_pt,bot_pt,glide_dir
        real(kind=DBL)                  :: z,new_length
        real(kind=DBL),dimension(3)     :: adj_length
        real(kind=DBL),dimension(2)     :: length

        ! For each of the two contacts...
        do j=1,2

            ! If contact is with the bottom..
            if (contacts(i,j,1) .lt. 1) then
                glide_dir(j,1:3)=0.d0
            else
                ! Determine the downward direction of the contacted fiber....
                top_pt(j,1:3)=fibers(int(contacts(i,j,1)),1:3)
                bot_pt(j,1:3)=fibers(int(contacts(i,j,1)),5:7)
                if (top_pt(j,3) .lt. bot_pt(j,3)) then
                    bot_pt(j,1:3)=top_pt(j,1:3)
                    top_pt(j,1:3)=fibers(int(contacts(i,j,1)),5:7)
                end if
                ! Get unit vector in downward contacted fiber direction, take z as slope...
                glide_dir(j,1:3)=(bot_pt(j,1:3)-top_pt(j,1:3))/fibers(int(contacts(i,j,1)),9)
            end if
            
            if (.true.) write(*,*) 'FIBER_GLIDE: glide_dir: ', glide_dir(j,1:3)
            !read(*,*)
            ! set z to positive for "along down direction"...
            z=0.d0-glide_dir(j,3)
            ! if z is below a threshold, no glide, stop further gliding...
            if ((z .lt. friction)) then
                glide_dir(j,1:3)=0.d0
                stop_glide(j)=.true.
            else
                if (debug) then
                    write(*,*) ' FIBER_GLIDE: Gliding! '
                    !read(*,*)
                end if
                stop_glide(j)=.false.
                ! if big z, set glide to x-y dir only of contacted fiber, the z fraction of a step.
                glide_dir(j,3)=0.d0
                glide_dir(j,1:3)=glide_dir(j,1:3)/norm2(glide_dir(j,1:3))*z*step_size
            end if
            if (.false.) then
                write(*,*) ' glide_dir: ',glide_dir(j,1:3)
                write(*,*) ' fibers(contacts(i,j,1),1:3), ...5:7): ',fibers(int(contacts(i,j,1)),1:3), &
                    fibers(int(contacts(i,j,1)),5:7)
                write(*,*) ' int(), contacts(i,j,1): ', int(contacts(i,j,1)), contacts(i,j,1)
                read (*,*)
            end if
        end do

        !!!THis should work for bottom touching, but is not quite right.
        !!! This is all unnecessary for both stop_glides=.true....

        one=2
        ! Determine which contact goes with which fiber end
        length(1)=norm2(contacts(i,1,2:4)-fibers(i,1:3))
        length(2)=norm2(contacts(i,2,2:4)-fibers(i,1:3))
        ! fiber end "one" (1:3) is closer to either contact 2 (one=2) or 1 (one=1)
        if (length(1) .lt. length(2)) one=1

        ! Set up new fiber position... WITH extra length...
        new_pos(:)=fibers(i,:)
        new_pos(1:3)=new_pos(1:3)+glide_dir(one,1:3)
        new_pos(5:7)=new_pos(5:7)+glide_dir(mod(one,2)+1,1:3)

        ! Compute new length, amount extended as vector along new fiber
        new_length=norm2(new_pos(1:3)-new_pos(5:7))
        adj_length(1:3)=(new_length-new_pos(9))*(new_pos(5:7)-new_pos(1:3))/new_length
        ! shorten each end by 0.5*extended vector....
        new_pos(5:7)=new_pos(5:7)-0.5*adj_length(1:3)
        new_pos(1:3)=new_pos(1:3)+0.5*adj_length(1:3)
        ! Check shortened new length...
        if (abs(norm2(new_pos(1:3)-new_pos(5:7))-new_pos(9)) .gt. 0.001) then
            write(*,*) ' GLIDE lenght ERROR: ',adj_length, new_pos(5:7)-new_pos(1:3), new_pos(9)
            read (*,*)
        end if

        if (.true.) then
            write(*,*) ' orig fiber: ',fibers(i,:)
            write(*,*) ' new_pos:    ',new_pos(:)
        end if


    end function




end module fibers_place

