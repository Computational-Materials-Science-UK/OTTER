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
    use bspline_module
    use bspline_kinds_module, only: wp, ip
    

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

        real(kind=DBL)                  ::  min_rad, max_rad, dir_deg, &
                                            big_box_vol,radius_range,length_range, &
                                            dir_rand, min_length, max_length, height, &
                                            max_height, friction
        real(kind=DBL),dimension(6)     ::  box_boundary,bigbox
        real(kind=DBL),dimension(3)     ::  box_range, shift
        integer                         ::  i,j,k, periodic_added, periodic_j, periodic_i, &
                                            max_fibers, go_glide, &
                                            scr_unit,dat_unit,rve_scr_unit,rve_nnc_unit, &
                                            rve,num_bad,excess_num,bending_case,union
        character(len=200)              ::  full_scr_name,full_nnc_name,full_name, &
                                            full_rves_batch_name,full_nncd_name, &
                                            saveas,exportstl
        character(len=50)               ::  date_time
        character(len=8)                ::  num_char
        logical                         ::  stopped, in_box, stop_slide
        logical,dimension(2)            ::  stop_glide
        
    

        !!!! Begin code...

        
        write (out_unit,*) ' Welcome to Fibers vb.1!'
        write (out_unit,*) ' '

        ! Get Global/General input, see otter_input
        call get_gen_input()

        ! Get fiber specific input...
        call get_fibers_input(min_rad, max_rad, min_length, max_length, friction)
        ! Adjust inputs...
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
        max_fibers=4
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

        ! contacts:  list of fiber contacts:
        !       1   - contacted fiber number, 
        !       2:4 - contact_pt3D (x,y,z), 
        !       5   - separation distance,
        !       6:8 - contact normal vector, NOT NORMALIZED
        !       9   - contact_case: 0 (side-to-side), 1 (contacting end to contacted side),
        !             2 (contacting side to contacted end), 3 (end-to-end)
        ! Each fiber has 2 (3 with rolling) plus those on top contacts..., 
        allocate(contacts(max_fibers,7,9))
        allocate(num_contacts(max_fibers))
        allocate(bent(max_fibers))
        allocate(fiber_path(max_fibers,41,3))
        allocate(control_points(max_fibers,5,3))
        allocate(bent_to_stright(max_fibers,41,1))
        bent=0

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
            !excess_num=1000000
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

                if (.false.) then
                    write (*,*) ' Enter fiber ',i,' by hand: x1, y1, dir_deg'
                    read(*,*) fibers(i,1), fibers(i,2), dir_deg
                    fibers(i,5)=fibers(i,1)+fibers(i,9)*cos(dir_deg/360.d0*2.d0*PI)
                    fibers(i,6)=fibers(i,2)+fibers(i,9)*sin(dir_deg/360.d0*2.d0*PI)
                end if
                !<<<<<<<<<<<<<<<<<<<<<<<<ADDED FOR TESTING>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                if (i.eq.1)then
                fibers(i,1)=1.d0
                fibers(i,2)=1.d0
                                           
                fibers(i,5)=fibers(i,1)
                fibers(i,6)=fibers(i,2)+fibers(i,9)
                end if
                if(i.eq.10)then
                    fibers(i,1)=10.d0
                fibers(i,2)=1.d0
                                           
                fibers(i,5)=fibers(i,1)
                fibers(i,6)=fibers(i,2)+fibers(i,9)
                end if
                if (i.eq.19)then
                    fibers(i,1)=25.d0
                    fibers(i,2)=1.d0
                                               
                    fibers(i,5)=fibers(i,1)
                    fibers(i,6)=fibers(i,2)+fibers(i,9)
                end if
                if (i.eq.28)then
                    fibers(i,1)=0.d0
                    fibers(i,2)=20.d0
                                               
                    fibers(i,6)=fibers(i,2)
                    fibers(i,5)=fibers(i,1)+fibers(i,9)
                end if

                

                
                !<<<<<<<<<<<<<<<<<<<<<<<<ADDED FOR TESTING>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                
                ! Drop fiber until contact, or hits bottom.
                stopped=.false.
                if (debug) write(out_unit,*) 'stop_glide: ',stop_glide(:)
                do while (.not. stopped)

                    ! Check for contacts
                    stop_slide=.false.
                    stop_glide(:)=.false.
                    call fiber_contact(i,max_rad,stop_slide,stop_glide)
                    if (debug) write(*,*) 'Re-entered program fibers_place from fiber_contact'

                    ! MOVE FIBER!
                    fibers(i,:)=fiber_move(i,stop_glide,stop_slide,friction,stopped)
                    
                end do ! dropping loop
                if (contacts(i,1,9).eq.4)bending_case=1
                if (contacts(i,1,1).gt.0.d0) then
                    call bending_test(i,bending_case)
                    bent(i)=1
                end if

                ! check if fiber is above box...
                height=max(fibers(i,3),fibers(i,7))
                if (height .gt. max_height) max_height=height
                if (min(fibers(i,3),fibers(i,7)) .gt. box_range(3)) num_bad=num_bad+1
                
                ! check fiber length...
                if (bent(i).eq.0)then
                    if (norm2(fibers(i,1:3)-fibers(i,5:7))-fibers(i,9) .gt. 0.001*max_length) then
                        write(*,*) ' FIBER LENGTH error! ',i
                        read(*,*)
                    end if
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
                                if (bent(i-periodic_added).eq.1)then
                                    control_points(i,:,:)=control_points(i-periodic_added,:,:)
                                    fiber_path(i,:,:)=fiber_path(i-periodic_added,:,:)
                                    control_points(i,:,1)=control_points(i,:,1)+(periodic_i-2)*box_range(1)
                                    fiber_path(i,:,1)=fiber_path(i,:,1)+(periodic_i-2)*box_range(1)
                                    control_points(i,:,2)=control_points(i,:,2)+(periodic_j-2)*box_range(2)
                                    fiber_path(i,:,2)=fiber_path(i,:,2)+(periodic_j-2)*box_range(2)
                                    bent(i)=1
                                end if
                            end if
                        end do
                    end do
                
                if (.false.) then
                    write(*,*) ' Do you want to stop? 1=yes'
                    read(*,*) j
                    if (j .eq. 1) max_fibers=i
                end if

                i=i+1
                
            end do ! fiber generation loop

            if (debug) write(out_unit,*) 'generating fibers complete'
            
            ! Check if sphere_num was too low...
            if (num_bad .le. excess_num) then
                write(out_unit,*) ' WARNING: excess_num not well selected - not a critical error.'
                !stop
            end if

            ! Write RVE to batch scr file...
            ! Write and union se(*,*)pheres...
            
            202 format('_sphere ',f0.3,',',f0.3,',',f0.3,' ',f0.3)
            203 format('_cone ',f0.3,',',f0.3,',',f0.3,' ',f0.3,' T ',f0.3,' A ',f0.3,',',f0.3,',',f0.3)
            201 format('_five ',f0.3,',',f0.3,',',f0.3,' ',f0.3,',',f0.3,',',f0.3,' '&
            ,f0.3,',',f0.3,',',f0.3,' ',f0.3,',',f0.3,',',f0.3,' ',f0.3,','&
            ,f0.3,',',f0.3,' ',f0.3)
            write(scr_unit,'(a)') '-osnap off'
            write(scr_unit,'(a)') '_zoom -300,-300,-300 500,500,500'
            union=0
            do j=1,i-1

                ! Check whether fiber has a least one end in box...
                in_box=.true.
                do k=1,3
                    if ( (min(fibers(j,k),fibers(j,k+4)) .gt. box_range(k)+max_rad) .or. &
                         (max(fibers(j,k),fibers(j,k+4)) .lt. 0.d0-max_rad) ) in_box=.false.
                end do

                ! If in the box, write to script...
                if (bent(j).eq.0)then
                    if (in_box) write(scr_unit,203) fibers(j,1:4),fibers(j,8),fibers(j,5:7)
                    if(debug)write(*,*)'[Otter_fibers 352] fiber number',j
                else
                    if (in_box) write(scr_unit,201) control_points(j,1,1:3),control_points(j,2,1:3),&
                    control_points(j,3,1:3),control_points(j,4,1:3),control_points(j,5,1:3), fibers(j,4)
                end if
                if (j.eq.50*union) then
                    write(scr_unit,'(a)') '_union all '
                    union=union+1
                end if
                 
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

        end do ! rves loop

        ! Clean up RVE scr...
        write(scr_unit,'(a)') 'FILEDIA 1'
        close(unit=scr_unit)
        

    end subroutine otter_fibers

    function fiber_rotate(i) result(delta_pos)
        implicit none

        ! Declare calling parameters:
        integer, intent(in)             :: i

        ! Output
        real(kind=DBL),dimension(9)     :: delta_pos

        !Declare local variables
        real(kind=DBL),dimension(3)     :: fiber_vec
        real(kind=DBL)                  :: theta,contact_ptR,r,f_perp,torque,i_rot, &
                                           dtheta, test_mid

        if(debug) write(out_unit,*) 'fiber_rotate SUBROUTINE: ',i,fibers(i,:)

        delta_pos(:)=0.d0
        
        ! Check for vertical fiber error...
        if (norm2(fibers(i,1:2)-fibers(i,5:6)) .le. 4.d0*step_size) then
            delta_pos(3)=0.d0 - 100000.d0
            delta_pos(7)=0.d0 - 100000.d0
            write(*,*) ' Dropping OUT! - LIKELY ERROR'
            read(*,*)
        else

            ! Find contact point on centerline:
            fiber_vec(:)=(fibers(i,5:7)-fibers(i,1:3))/fibers(i,9)
            contact_ptR=dot_product((contacts(i,1,2:4)-fibers(i,1:3)),fiber_vec(:))

            ! Check for rocking on mid-point...
            test_mid=contact_ptR-0.5*fibers(i,9)
            if (abs(test_mid) .le. step_size) then
                fibers(i,1:3)=fibers(i,1:3)+sign(1.d0,-1.d0*test_mid)*step_size*fiber_vec(:)
                fibers(i,5:7)=fibers(i,5:7)+sign(1.d0,-1.d0*test_mid)*step_size*fiber_vec(:)
                contact_ptR=dot_product((contacts(i,1,2:4)-fibers(i,1:3)),fiber_vec(:))
                !write(*,*) ' FIBER_ROTATE: shifting fiber off mid-point... '
            end if

            ! Find r vector for torque calculation...
            r=contact_ptR-0.5*fibers(i,9)
            ! Compute force F_perp=F_g*cos(asin(z/1))
            f_perp=0.d0-cos(asin(fiber_vec(3)))
            ! torque is just r*f_perp, maintaining signs for sense...
            torque=r*f_perp

            ! Get rotation moment of inertia, strictly this is I/mass...
            ! This is i_rot from contact to end1 PLUS i_rot from contact to end2, both are solid cylinder i_rot's, with radius fixed to average radius....  **KLUDGE**
            i_rot=0.5*((fibers(i,4)+fibers(i,8))/2.d0)**2+( (fibers(i,9)-contact_ptR)**2 + (contact_ptR)**2 )/3.d0
            ! Get dtheta and change in z!
            dtheta=torque/i_rot
            delta_pos(3)=dtheta*contact_ptR
            delta_pos(7)=dtheta*( contact_ptR - fibers(i,9) )

            ! Switch rotation direciton if 5:7 is long end and goes up...
            if ((r .lt. 0.d0) .and. (delta_pos(7) .gt. 0.d0)) delta_pos(:)=0.d0-delta_pos(:)

            if (.false.) then
                write(*,*) ' ROTATE ................'
                write(*,*) ' contact_ptR, f_perp, torque: ',contact_ptR, f_perp, torque
                write(*,*) ' i_rot, dtheta, delta_pos: ',i_rot, dtheta, delta_pos(3), delta_pos(7)
                read(*,*)
            end if

        end if

    end function fiber_rotate


    subroutine fiber_contact(i,max_rad,stop_slide,stop_glide)
        implicit none 

        !Declare calling parameters
        integer, intent(in) :: i
        real(kind=DBL),intent(in) :: max_rad
        logical,intent(inout) :: stop_slide
        logical,dimension(2),intent(inout) :: stop_glide

        !Declare local variables
        real(kind=DBL), dimension(3)    :: contact_vec,contact_pt3D,contact_norm
        real(kind=DBL), dimension(2)    :: radius,near_ptP,sign_contact, near_endP
        real(kind=DBL), dimension(2,3)  :: near_pt3D,fiber_vec,face_vec, &
                                           center, near_end3D
        real(kind=DBL)                  :: dist,bottom, dist_end, sign_fiber
        integer                         :: j,k,contact_case,l,free_end,stopped_end
        logical                         :: contact,alt_contact=.false.
        !real(kind=DBL), dimension(max_fibers,1) :: bent
        !bent=0.d0
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
                if (bent(j).eq.1) then
                   
                    call bent_contact(i,j,contact,contact_pt3D)
                    !write(*,*)'Here 1'
                else
                    if (debug) write (*,*) ' Check Fine contact: ',i,j
                    contact=.false.
                    !write(*,*)'Here 2'
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
                                contact_norm(1:3)=sign_fiber*fiber_vec(2,:)
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
                            
                            l=j
                            if (k .eq. 2) l=i
                            call segments_dist_3d(center(k,:),center(k,:)+face_vec(k,:),fibers(l,1:3), &
                                                fibers(l,5:7),dist_end,near_end3D(1,:), &
                                                near_end3D(2,:), near_endP(1), near_endP(2) )
                            if (dist_end .lt. radius(mod(k,2)+1)-min_olp) then
                                contact_pt3D(1:3)=near_end3D(1,:)
                                if (k .eq. 2) then 
                                    contact_norm(1:3)=sign_fiber*fiber_vec(k,:)
                                else
                                    contact_norm(1:3)=near_end3D(1,:)-near_end3D(2,:)
                                end if
                                contact=.true. 
                            end if

                        case default

                            !!!! ORIGINAL ROUNDED END
                            if (dist .le. (radius(1)+radius(2)-min_olp)) then 
                                contact=.true.
                                contact_pt3D(1:3)=near_pt3D(1,1:3)+0.5*(near_pt3D(2,1:3)-near_pt3D(1,1:3))
                                contact_norm(1:3)=near_pt3D(1,:)-near_pt3D(2,:)
                            end if
                            !!!! END ORIGINAL ROUNDED END

                        end select

                    end if
                end if  
            else
                
                

            end if ! end fine check
            !write(*,*)'Here 612'
            if (contact) then
                !write(*,*) ' Record fine contact! ',i,j
                !write(*,*) ' Fiber 1: ',fibers(i,1:3),fibers(i,5:7)
                !write(*,*) ' Fiber 2: ',fibers(j,1:3),fibers(j,5:7)
                !write(*,*) ' dist, near_ptP(1),(2): ',dist, near_ptP(1), near_ptP(2)
                !write(*,*) ' near_pt3D(1),(2): ',near_pt3D(1,:), near_pt3D(2,:) 
                num_contacts(i)=num_contacts(i)+1
                contacts(i,num_contacts(i),1)=j
                contacts(i,num_contacts(i),2:4)=contact_pt3D(1:3)
                contacts(i,num_contacts(i),5)=dist
                contacts(i,num_contacts(i),6:8)=contact_norm(1:3)/norm2(contact_norm(1:3))
                contacts(i,num_contacts(i),9)=contact_case
                !write(*,*)'Here 625'
                ! If end of falling fiber is contact, set stop_slide and stop_glide, if appropriate
                if ((contact_case .eq. 1) .or. (contact_case .eq. 3)) then
                    ! Check if falling fiber lower end is the contact end... if so, stop_slide
                    stopped_end=4*int(near_ptP(1))+3 ! 0->3, 1->7
                    free_end=7-4*int(near_ptP(1)) ! 0->7, 1->3
                    if (fibers(i,stopped_end) .lt. fibers(i,free_end)) stop_slide=.true.
                    if (near_ptP(1) .eq. 0.d0) then
                        stop_glide(1)=.true.
                    else
                        stop_glide(2)=.true.
                    end if
                end if
            
                if (.false.) then
                    write(*,*) ' Record fine contact!, Contact Case ',i,j,contact_case
                    write(*,*) ' Contacts: ',contacts(i,num_contacts(i),:)
                end if
            end if

        end do

        bottom=0.d0-0.25*box_length(3)
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TESTING>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        bottom=1.d0!5.d-1
        !!!!!!!!  SHOULD ONLY BE FOR TEST_GLIDE....
        !bottom=0
        do k=3,7,4
            if (fibers(i,k) .le. bottom) then
                num_contacts(i)=num_contacts(i)+1
                contacts(i,num_contacts(i),1)=0
                contacts(i,num_contacts(i),2:4)=fibers(i,k-2:k)
                contacts(i,num_contacts(i),5)=0.d0
                contacts(i,num_contacts(i),6:8)=[0.d0,0.d0,1.d0]
                contacts(i,num_contacts(i),9)=4
                stop_slide=.true.
                stop_glide(int(k/4)+1)=.true.
                write(*,*) ' One end is touching the bottom...', k
                !read(*,*)
            end if
        end do

    end subroutine

    function fiber_glide(i, friction, stop_glide) result(delta_pos)
        implicit none

        ! Declare calling parameters:
        integer, intent(in)             :: i
        real(kind=DBL),intent(in)       :: friction

        ! Output
        real(kind=DBL),dimension(9)     :: delta_pos
        logical,dimension(2)            :: stop_glide

        ! Local variables...
        integer                         :: j,one
        real(kind=DBL),dimension(2,3)   :: top_pt,bot_pt,glide_dir
        real(kind=DBL)                  :: residual_force_frac,contact_ptR
        real(kind=DBL),dimension(2)     :: length,dist
        real(kind=DBL),dimension(3)     :: fiber_vec

        if (.false.) write(*,*) ' Entering FIBER_GLIDE...'

        ! Get fiber_vec... this is normalized...
        fiber_vec(:)=(fibers(i,5:7)-fibers(i,1:3))/fibers(i,9)

        ! For each of the two contacts...
        do j=1,2
            ! If contact is with the bottom..
            if (contacts(i,j,9) .eq. 4) then
                glide_dir(j,1:3)=0.d0
                dist(j)=0.5*fibers(i,9)
                stop_glide(j)=.true.
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

                ! Find contact point dist to midpoint:
                contact_ptR=dot_product((contacts(i,j,2:4)-fibers(i,1:3)),fiber_vec(:))
                dist(j)=abs(contact_ptR-0.5*fibers(i,9))

                ! Add friction...
                if (glide_dir(j,3) .ne. 0.d0) then
                    residual_force_frac=1.d0-friction/abs(glide_dir(j,3))
                else
                    residual_force_frac=0.d0
                end if
                glide_dir(j,3)=residual_force_frac*glide_dir(j,3)
                if (glide_dir(j,3) .ge. 0.d0) then
                    glide_dir(j,:)=0.d0
                    stop_glide(j)=.true.
                end if
            end if
        end do

        ! With both dist(j)'s calculated, get the real glides...
        do j=1,2
            
            if (.not. stop_glide(j)) then
                ! Weight by load at contact
                glide_dir(j,1:3)=glide_dir(j,1:3)*dist(mod(j,2)+1)/(dist(1)+dist(2))
            else
                glide_dir(j,1:3)=0.d0
            end if
            
            if (.false.) then
                write(*,*) ' FIBER_GLIDE: j,glide_dir: ',j,glide_dir(j,1:3)
                !write(*,*) ' fibers(contacts(i,j,1),1:3), ...5:7): ',fibers(int(contacts(i,j,1)),1:3), &
                !    fibers(int(contacts(i,j,1)),5:7)
                !write(*,*) ' int(), contacts(i,j,1): ', int(contacts(i,j,1)), contacts(i,j,1)
                !read (*,*)
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
        delta_pos(:)=0.d0
        delta_pos(1:3)=delta_pos(1:3)+glide_dir(one,1:3)
        delta_pos(5:7)=delta_pos(5:7)+glide_dir(mod(one,2)+1,1:3)

        if (.false.) then
            !write(*,*) ' orig fiber: ',fibers(i,:)
            write(*,*) ' FIBER_GLIDE: delta_pos:    ',delta_pos(:), stop_glide(:)
        end if


    end function


    function fiber_slide(i,stop_slide,friction) result(new_pos)
        implicit none

        ! Declare calling parameters:
        integer, intent(in)             :: i
        real(kind=DBL),intent(in)       :: friction

        ! Output
        real(kind=DBL),dimension(9)     :: new_pos
        logical, intent(inout)          :: stop_slide

        ! Local variables
        real(kind=DBL),dimension(3)     :: top_pt,bot_pt,fiber_vec
        real(kind=DBL)                  :: fiber_angle_ratio, residual_slide_force
        integer                         :: j

        !!!! *** DOES NOT SLIDE OFF ENDS!!!!!

        if (.false.) write(*,*) ' Entering FIBER_SLIDE...'
        ! Default is no motion!...
        new_pos(:)=0.d0

        if (.not. stop_slide) then
            !! Get angle of fiber
            ! Get fiber vector in DOWNWARD direction, compute sin(angle)
            top_pt(1:3)=fibers(i,1:3)
            bot_pt(1:3)=fibers(i,5:7)
            if (top_pt(3) .lt. bot_pt(3)) then
                bot_pt(1:3)=top_pt(1:3)
                top_pt(1:3)=fibers(i,5:7)
            end if
            fiber_vec(:)=bot_pt(:)-top_pt(:)
            ! fiber_angle_ratio=rise (opposite) / length (hypo) = sin(angle)
            fiber_angle_ratio=( top_pt(3)-bot_pt(3) )/norm2(fiber_vec(:))
            ! Slide force is F*sin(angle) - F*sin(friction_angle), for F=1, friction=asin(friction_angle)
            residual_slide_force=fiber_angle_ratio-friction
        
            !! If greater than zero, then move in downward direction by step_size...
            if (residual_slide_force .gt. 0.d0) then
                ! get fiber direction downwards...
                fiber_vec(:)=fiber_vec(:)/norm2(fiber_vec(:))
                new_pos(1:3)=residual_slide_force*fiber_vec(:)
                new_pos(5:7)=residual_slide_force*fiber_vec(:)
            else
                if (.false.) write(*,*) ' FIBER_SLIDE: Too shallow, stopping SLIDE.'
                stop_slide=.true.
            end if
            if (.false.) write(*,*) ' FIBER_SLIdE: new_pos: ',new_pos(:)
        end if

    end function fiber_slide

    function fiber_move(i,stop_glide,stop_slide,friction,stopped) result(new_pos)
        implicit None

        ! Inputs
        integer                     :: i
        logical                     :: stop_slide,stopped
        logical,dimension(2)        :: stop_glide
        real(kind=DBL)              :: friction

        ! Outputs
        real(kind=DBL),dimension(9) :: new_pos

        ! Local variables
        real(kind=DBL)              :: max_move, alt_move, new_length
        real(kind=DBL),dimension(3) :: adj_length

        !! MECHANICAL MODEL: We will effectively time step motion by assuming that real forces pull each fiber end. Take m*g=1 and t^2=2*m to set the scale, which is t=sqrt(2/g). With these assumptions, the motion due to gravity of the fiber per time step is 1 unit.  This is then normalized by the "step_size".

        !! For slides, the force is F_g*sin(theta) along the fiber, with theta with respect to XY plane.  Therefore, in one time step t=sqrt(2/g) we have motion of 0.5*m*g*sin(theta)/m*[sqrt(2/g)]^2=sin(theta)... Both ends move sin(theta) along the fiber.

        !! Glides are like slides, but each end slides independently. For an end i, and the contact nearest that end, ci (at point ciP from the first end), the distance to the center of gravity of the fiber is d_ci. The other end is j, cj, cjP, and d_cj.  The F_z on contact i is therefore F_zi=F_g*d_ci/(d_ci+d_cj).  The glide force for that end is therefore F_zi*sin(theta), for theta the angle fo the contacted fiber. The motion of that point is d_ci/(d_ci+d_cj)*sin(theta)... Use the new ci and cj with the OLD ciP and cjP to compute the NEW end points i and j.

        !! For rotations, we have the moment M=(contact_point - cog) CROSS z_hat.  For t=sqrt(2/g) we have that M/I=d(theta) and dzi=d(theta)*Li, for Li the distance of end i from the contact point.  Here, I=0.5*r^2+1/3*(L1^2+L2^2)

        !! This gives us a dzi from dropping (=1, when applicable) and rotation (where applicable), plus a d(xyz)i from sliding and gliding.  Summing all of these gives a dri, which we then normalize so that the largest equals the step size.

        !! apply to fiber, renormalize length

        if (.false.) write(*,*) ' fiber_move: num_contacts: ',num_contacts(i)

        ! Get the raw move...
        new_pos(:)=get_move(i,stop_glide,stop_slide,friction,stopped)
        if (.false.) write (*,*) new_pos(:)

        ! Find end that moves the most, and normalize that motion to step_size
        max_move=norm2(new_pos(1:3))
        alt_move=norm2(new_pos(5:7))
        if (alt_move .gt. max_move) max_move=alt_move

        if (max_move .ne. 0.d0) then
            new_pos(1:3)=new_pos(1:3)/max_move*step_size
            new_pos(5:7)=new_pos(5:7)/max_move*step_size

            ! With 1:3,5:7 the change in position, and everything else zero, we can add fibers(i,:)
            new_pos(:)=new_pos(:)+fibers(i,:)
            if (.false.) write (*,*) ' FIBER_MOVE: new_pos: ',new_pos(:)

            ! Normalize length
            ! Compute new length, amount extended as vector along new fiber
            new_length=norm2(new_pos(1:3)-new_pos(5:7))
            adj_length(1:3)=(new_length-new_pos(9))*(new_pos(5:7)-new_pos(1:3))/new_length
            ! shorten each end by 0.5*extended vector....
            new_pos(5:7)=new_pos(5:7)-0.5*adj_length(1:3)
            new_pos(1:3)=new_pos(1:3)+0.5*adj_length(1:3)
            ! Check shortened new length...
            if (abs(norm2(new_pos(1:3)-new_pos(5:7))-new_pos(9)) .gt. 0.001) then
                write(*,*) ' MOVE lenght ERROR: ',adj_length, new_pos(5:7)-new_pos(1:3), new_pos(9)
                read (*,*)
            end if
        else
            new_pos(:)=fibers(i,:)
        end if

    end function fiber_move

    recursive function get_move(i,stop_glide,stop_slide,friction,stopped) result(delta_pos)

        ! Inputs
        integer                     :: i
        logical                     :: stop_slide,stopped
        logical,dimension(2)        :: stop_glide
        real(kind=DBL)              :: friction

        ! Outputs
        real(kind=DBL),dimension(9) :: delta_pos

        ! Local variables
        logical                     :: overhang
        integer                     :: j,k,overhang_contact
        real(kind=DBL)              :: sum_top,sum_bot
        real(kind=DBL),dimension(3) :: moments,mid_point
        real(kind=DBL),dimension(9) :: temp
        
        delta_pos(:)=0.d0

        select case (num_contacts(i))
        case (0)
            ! Drop fiber by step_size
            delta_pos(3)=delta_pos(3)-1
            delta_pos(7)=delta_pos(7)-1
        case (1)
            if (.not. stop_slide) delta_pos(:)=delta_pos(:)+fiber_slide(i,stop_slide,friction)
            delta_pos(:)=delta_pos(:)+fiber_rotate(i)
        case (2)
            ! Check overhang...
            overhang=check_overhang(i,overhang_contact)
            if (.false.) write(*,*) ' Overhang check is: ',overhang
            if (overhang) then
                num_contacts(i)=1
                contacts(i,1,:)=contacts(i,overhang_contact,:)
                delta_pos(:)=delta_pos(:)+get_move(i,stop_glide,stop_slide,friction,stopped)
            else ! If not overhung, then SLIDE and GLIDE!

                ! Call fiber slide with double the friction!
                if (.not. stop_slide) delta_pos(:)=delta_pos(:)+fiber_slide(i,stop_slide,friction*2.d0)

                ! IF at least one stop_glide is false, then start gliding...
                if (.not. (stop_glide(1) .and. stop_glide(2))) then
                    delta_pos(:)=delta_pos(:)+fiber_glide(i,friction,stop_glide)

                    ! IF fiber_glide sets both stop_glides to true, then can it still slide?
                    if (stop_glide(1) .and. stop_glide(2)) then
                        if (stop_slide) then
                            stopped=.true.
                            if (.false.) write(out_unit,*) ' Stopped with two contacts, both below friction, and stop_slide. '
                        end if
                    else  ! else, keep gliding...
                        if (.false.) write(out_unit,*) ' Keep gliding...'
                    end if
                else ! Check if sliding is allowed...
                    if (stop_slide) then
                        if (.false.) write(out_unit,*) ' stopped with two stop_glides and stop_slide.'
                        !read(*,*)
                        stopped=.true.
                    end if
                end if
            end if

            if (.false.) write(*,*) ' Done with two contact move.'

        case default
            !! Find mid-point for moments...
            mid_point(:)=(fibers(i,1:3)+fibers(i,5:7))/2

            !! Sort contacts bottom to top...
            do j=1,num_contacts(i)-1
                do k=j+1,num_contacts(i)
                    if (contacts(i,j,4) .gt. contacts(i,k,4)) then
                        temp(:)=contacts(i,j,:)
                        contacts(i,j,:)=contacts(i,k,:)
                        contacts(i,k,:)=temp(:)
                    end if
                end do
            end do

            !! Compute and "normalize" moments...
            do j=1,num_contacts(i)
                moments(j)=cross2(( contacts(i,j,2:3)-mid_point(1:2) ),contacts(i,j,6:7) )
                if (moments(j) .gt. 0.d0) then
                    moments(j)=1.d0
                else
                    if (moments(j) .lt. 0.d0) moments(j)=-1.d0
                end if
            end do

            !! Check if bottom contact is touching the bottom...
            if (contacts(i,1,9) .eq. 4.d0) then
                if (moments(2)+moments(3) .eq. 0.d0) then
                    stopped=.true.
                else
                    ! find contact with shallowest slope, therefore steepest normal...
                    if (contacts(i,2,8) .lt. contacts(i,3,8)) contacts(i,2,:)=contacts(i,3,:)
                    num_contacts(i)=num_contacts(i)-1
                    stop_glide(1)=.true.
                    delta_pos(:)=delta_pos(:)+get_move(i,stop_glide,stop_slide,friction,stopped)
                end if
            else
                sum_top=0.d0
                sum_bot=0.d0
                do j=1,num_contacts(i)
                    if (contacts(i,j,4) .gt. mid_point(3)) sum_top=sum_top+moments(j)
                    if (contacts(i,j,4) .lt. mid_point(3)) sum_bot=sum_bot+moments(j)
                end do
                if ((sum_top .eq. 0.d0) .or. (sum_bot .eq. 0.d0)) then
                    stopped=.true.
                else
                    do j=1,num_contacts(i)-1
                        ! Put least steep normal, highest slope, in last position to be ignored
                        if (contacts(i,num_contacts(i),8) .lt. contacts(i,j,8)) then
                            temp(:)=contacts(i,num_contacts(i),:)
                            contacts(i,num_contacts(i),:)=contacts(i,j,:)
                            contacts(i,j,:)=temp(:)
                        end if
                    end do
                    num_contacts(i)=num_contacts(i)-1
                    delta_pos(:)=delta_pos(:)+get_move(i,stop_glide,stop_slide,friction,stopped)
                end if
            end if

            if (.false.) write(*,*) ' Done with three contact move.'

            !if (.true.) write(out_unit,*) ' Stopped w/ 3 or more contacts...',contacts(i,1:3,5)
            !stopped=.true.
        end select

    
    end function get_move

    function check_overhang(i,overhang_contact) result(overhang)
        implicit None

        ! Inputs
        integer,intent(in)              :: i
        integer,intent(inout)           :: overhang_contact

        ! Outputs
        logical                         :: overhang

        ! Local variables
        integer                         :: end1_close,end2_close,j,k
        real(kind=DBL),dimension(3)     :: mid_point
        real(kind=DBL),dimension(9)     :: temp,dist_end1,dist_end2
        real(kind=DBL)                  :: temp1

        end1_close=0
        end2_close=0
        overhang=.false.
        overhang_contact=0

        ! Compute fiber midpoint
        mid_point(:)=(fibers(i,1:3)+fibers(i,5:7))/2
        do j=1,num_contacts(i)
            dist_end2(j)=norm2(contacts(i,j,2:4)-fibers(i,5:7))
            dist_end1(j)=norm2(contacts(i,j,2:4)-fibers(i,1:3))
            if (dist_end1(j) .lt. dist_end2(j)) then
                end1_close=end1_close+1
            else
                if (dist_end1(j) .gt. dist_end2(j)) then
                    end2_close=end2_close+1
                end if
            end if
        end do
        if ((end1_close .eq. 0) .or. (end2_close .eq. 0)) then
            overhang=.true.
            do j=1,num_contacts(i)-1
                do k=j+1,num_contacts(i)
                    if (dist_end1(j) .gt. dist_end1(k)) then
                        temp(:)=contacts(i,j,:)
                        contacts(i,j,:)=contacts(i,k,:)
                        contacts(i,k,:)=temp(:)
                        temp1=dist_end1(j)
                        dist_end1(j)=dist_end1(k)
                        dist_end1(k)=temp1
                        temp1=dist_end2(j)
                        dist_end2(j)=dist_end2(k)
                        dist_end2(k)=temp1
                    end if
                end do
            end do
            if (end2_close .eq. 0) then
                overhang_contact=num_contacts(i)
            else
                overhang_contact=1
            end if
        end if
        
    end function check_overhang
subroutine bent_contact(i,j,contact,contact_pt3D)
    integer:: max_fibers
    
    real(kind=8),dimension(3):: pt,endpoint1,endpoint2,AB,BE,AE,AB_cross_AE,AB_cross_BE
    real(kind=DBL),dimension(3)::p,temp_p,contact_pt3D,temp_pt
    real(kind=DBl),dimension(41)::slope_PE
    real(kind=8):: AB_dot_BE,AB_dot_AE,mag_AB,t1,t2,t3,temp_dist,temp_slope,AB_dot_AB 
    integer(ip)::i,n,k
    logical                         :: contact
    
     
    t=0
    temp_slope=0
    temp_dist=10.d0
   
    endpoint1(1:3)=fibers(i,1:3)
    endpoint2(1:3)=fibers(i,5:7)
    ! endpoint1(1)=0.d0
    ! endpoint1(2)=0.d0
    ! endpoint1(3)=0.d0
    ! endpoint2(1)=10.d0
    ! endpoint2(2)=0.d0
    ! endpoint2(3)=5.d0
    AB=endpoint2-endpoint1
    AB_dot_AB=dot_product(AB,AB)
    !write(*,*)'End Point 1 = ',endpoint1
    !write(*,*)'End Point 2 = ',endpoint2
    ! write(*,*)'AB = ',AB
    do n=1,21
        ! write(*,*)'j = ',j
        ! write(*,*)'n = ',n
    
        pt(1)=fiber_path(j,n,1)
        pt(2)=fiber_path(j,n,2)
        pt(3)=fiber_path(j,n,3)
        
        
        !write(*,*)'pt',pt
        

        
        BE=pt-endpoint2
        AE=pt-endpoint1
        
        ! write(*,*)'BE = ',BE
        ! write(*,*)'AE = ',AE
        AB_dot_BE=DOT_PRODUCT(AB,BE)
        AB_dot_AE=DOT_PRODUCT(AB,AE)
        if (AB_dot_AE.lt.0.d0) then
            bent_to_stright(k,n,1)=sqrt(AE(1)**2.d0+AE(2)**2.d0+AE(3)**2d0)
           
        else
        if (AB_dot_BE.gt.0.d0)then
            bent_to_stright(k,n,1)=sqrt(BE(1)**2.d0+BE(2)**2.d0+BE(3)**2d0)
            
        else
            AB_cross_AE(1)=AB(2)*AE(3)-(AB(3)*AE(2))
            AB_cross_AE(2)=AB(3)*AE(1)-(AB(1)*AE(3))
            AB_cross_AE(3)=AB(1)*AE(2)-(AB(2)*AE(1))
            !write(*,*)'AB X AE = ',AB_cross_AE
            mag_AB=sqrt(AB(1)**2.d0+AB(2)**2.d0+AB(3)**2.d0)
            !write(*,*)'mag.AB = ',mag_AB
            bent_to_stright(k,n,1)=sqrt((AB_cross_AE(1)**2.d0+AB_cross_AE(2)**2.d0+AB_cross_AE(3)**2.d0))/mag_AB
            
            t=AB_dot_AE/AB_dot_AB
            ! write(*,*)'T = ',t
            ! write(*,*)'P = ',p
            P=AB*t+endpoint1
            slope_PE(n)=(P(3)-pt(3))/sqrt((P(1)-pt(1))**2+(P(2)-pt(2))**2)
            

        end if
        end if
        !write(*,*)'[Bent Contact] fiber(K) =',k
        !write(*,*)'[Bent Contact] Point(n) =',n
        !write(*,*)'AB . BE = ',AB_dot_BE
        !write(*,*)'AB . AE = ',AB_dot_AE
        
        if(n.gt.1)then
            if(slope_PE(n).gt.temp_slope)then
                temp_slope=slope_PE(n)
                temp_p=P
                temp_pt=pt
                temp_dist=bent_to_stright(k,n,1)
            end if    
        endif
        
    end do
    if (temp_dist.lt.(fibers(i,4)+fibers(j,4)))then
        contact=.true.
        contact_pt3D(1:3)=temp_pt(1:3)+0.5*(temp_p(1:3)-temp_pt(1:3))
    else
        contact=.false.
    end if

        !write(*,*)'Slope = ',temp_slope
        !write(*,*)'Distance = ',temp_dist
        !write(*,*)'point ',temp_p
        !write(*,*)'Fiber point =',temp_pt
        !if(contact)read(*,*)
        
        
    end subroutine bent_contact

    subroutine bending_test(i,bending_case)
        implicit none
        integer, intent(in)                 :: i
        integer                             :: bending_case,num_pt
        real(kind=DBL)                      :: bottom
        real(kind=DBL),dimension(4)         :: p1,p2,endpoint1,endpoint2,cog
        real(wp),dimension(5)               :: dist
        
         
        
        num_pt=5
    !  bottom=0.d0-0.25*box_length(3)
    !  bottom=1.d0
     endpoint1=1.d0
     endpoint2=1.d0
     endpoint1(1:3)=fibers(i,1:3)
     endpoint2(1:3)=fibers(i,5:7)
     p1=1.d0
     p1(1:3)=endpoint1(1:3)+((endpoint2(1:3)-endpoint1(1:3))/4.d0)
     p2=1.d0
     p2(1:3)=endpoint1(1:3)+((endpoint2(1:3)-endpoint1(1:3))/(4.0/3.0))
     cog=1.d0
     cog(1:3)=endpoint1(1:3)+((endpoint2(1:3)-endpoint1(1:3))/2.d0)
     dist=0.d0

     if(num_contacts(i).eq.2)then
         if (norm2(contacts(i,1,2:4)-endpoint1(1:3)).lt.norm2(contacts(i,2,2:4)-endpoint1(1:3))) then        
              dist(1)=norm2(contacts(i,1,2:4)-endpoint1(1:3))
         else
              dist(1)=norm2(contacts(i,2,2:4)-endpoint1(1:3))
         end if
         if (norm2(contacts(i,1,2:4)-p1(1:3)).lt.norm2(contacts(i,2,2:4)-p1(1:3))) then        
            dist(2)=norm2(contacts(i,1,2:4)-p1(1:3))
       else
            dist(2)=norm2(contacts(i,2,2:4)-p1(1:3))
       end if
         if (norm2(cog(1:3)-contacts(i,1,2:4)).lt.norm2(cog(1:3)-contacts(i,2,2:4))) then
              dist(3)=norm2(cog(1:3)-contacts(i,1,2:4))
         else
              dist(3)=norm2(cog(1:3)-contacts(i,2,2:4))
         end if
         if (norm2(contacts(i,1,2:4)-p2(1:3)).lt.norm2(contacts(i,2,2:4)-p2(1:3))) then        
            dist(4)=norm2(contacts(i,1,2:4)-p2(1:3))
       else
            dist(4)=norm2(contacts(i,2,2:4)-p2(1:3))
       end if
         if(norm2(endpoint2(1:3)-p1(1:3)).lt.norm2(endpoint2(1:3)-p2(1:3))) then
                 dist(5)=norm2(endpoint2(1:3)-p1(1:3))
         else
                 dist(5)=norm2(endpoint2(1:3)-p2(1:3))
         end if
         !bent(i,1)=1
     end if
     if(num_contacts(i).eq.1)then 
             dist(1)=norm2(contacts(i,1,2:4)-endpoint1(1:3))
             dist(2)=norm2(contacts(i,1,2:4)-p1(1:3))
             dist(4)=norm2(contacts(i,1,2:4)-p2(1:3))
             dist(5)=norm2(contacts(i,1,2:4)-endpoint2(1:3))
             !bent(i,1)=1
     end if
     if(num_contacts(i).eq.3)then
        if (norm2(contacts(i,1,2:4)-endpoint1(1:3)).lt.norm2(contacts(i,2,2:4)-endpoint1(1:3))) then 
            if (norm2(contacts(i,1,2:4)-endpoint1(1:3)).lt.norm2(contacts(i,3,2:4)-endpoint1(1:3))) then
                dist(1)=norm2(contacts(i,1,2:4)-endpoint1(1:3))
            else
                dist(1)=norm2(contacts(i,3,2:4)-endpoint1(1:3))
            end if
        else
            if (norm2(contacts(i,2,2:4)-endpoint1(1:3)).lt.norm2(contacts(i,3,2:4)-endpoint1(1:3))) then
                dist(1)=norm2(contacts(i,2,2:4)-endpoint1(1:3))
            else
                dist(1)=norm2(contacts(i,3,2:4)-endpoint1(1:3))
            end if
        end if

        if (norm2(contacts(i,1,2:4)-p1(1:3)).lt.norm2(contacts(i,2,2:4)-p1(1:3))) then 
            if (norm2(contacts(i,1,2:4)-p1(1:3)).lt.norm2(contacts(i,3,2:4)-p1(1:3))) then
                dist(2)=norm2(contacts(i,1,2:4)-p1(1:3))
            else
                dist(2)=norm2(contacts(i,3,2:4)-p1(1:3))
            end if
        else
            if (norm2(contacts(i,2,2:4)-p1(1:3)).lt.norm2(contacts(i,3,2:4)-p1(1:3))) then
                dist(2)=norm2(contacts(i,2,2:4)-p1(1:3))
            else
                dist(2)=norm2(contacts(i,3,2:4)-p1(1:3))
            end if
        end if

        if (norm2(contacts(i,1,2:4)-cog(1:3)).lt.norm2(contacts(i,2,2:4)-cog(1:3))) then 
            if (norm2(contacts(i,1,2:4)-cog(1:3)).lt.norm2(contacts(i,3,2:4)-cog(1:3))) then
                dist(3)=norm2(contacts(i,1,2:4)-cog(1:3))
            else
                dist(3)=norm2(contacts(i,3,2:4)-cog(1:3))
            end if
        else
            if (norm2(contacts(i,2,2:4)-cog(1:3)).lt.norm2(contacts(i,3,2:4)-cog(1:3))) then
                dist(3)=norm2(contacts(i,2,2:4)-cog(1:3))
            else
                dist(3)=norm2(contacts(i,3,2:4)-cog(1:3))
            end if
        end if

        if (norm2(contacts(i,1,2:4)-p2(1:3)).lt.norm2(contacts(i,2,2:4)-p2(1:3))) then 
            if (norm2(contacts(i,1,2:4)-p2(1:3)).lt.norm2(contacts(i,3,2:4)-p2(1:3))) then
                dist(4)=norm2(contacts(i,1,2:4)-p2(1:3))
            else
                dist(4)=norm2(contacts(i,3,2:4)-p2(1:3))
            end if
        else
            if (norm2(contacts(i,2,2:4)-p2(1:3)).lt.norm2(contacts(i,3,2:4)-p2(1:3))) then
                dist(4)=norm2(contacts(i,2,2:4)-p2(1:3))
            else
                dist(4)=norm2(contacts(i,3,2:4)-p2(1:3))
            end if
        end if

        if (norm2(contacts(i,1,2:4)-endpoint2(1:3)).lt.norm2(contacts(i,2,2:4)-endpoint2(1:3))) then 
            if (norm2(contacts(i,1,2:4)-endpoint2(1:3)).lt.norm2(contacts(i,3,2:4)-endpoint2(1:3))) then
                dist(5)=norm2(contacts(i,1,2:4)-endpoint2(1:3))
            else
                dist(5)=norm2(contacts(i,3,2:4)-endpoint2(1:3))
            end if
        else
            if (norm2(contacts(i,2,2:4)-endpoint2(1:3)).lt.norm2(contacts(i,3,2:4)-endpoint2(1:3))) then
                dist(5)=norm2(contacts(i,2,2:4)-endpoint2(1:3))
            else
                dist(5)=norm2(contacts(i,3,2:4)-endpoint2(1:3))
            end if
        end if
    end if
     

     !call fiber_mapping(i,dist,endpoint1,endpoint2,cog,p1,p2)
     
     call fiber_mapping_alt(i,dist,endpoint1,endpoint2,cog,p1,p2)

        !Case 1: Touching the ground and one fiber
        if (bending_case.eq.1)then
            


        ! write(*,*)'[Bending Test] i =',i
        ! write(*,*)'[Bending Test] num_contact =',num_contacts(i)
        ! write(*,*)'[Bending Test] contact pt =',contacts(i,1,1)
        ! write(*,*)'[Bending Test] contact pt =',contacts(i,1,2:4)
        ! write(*,*)'[Bending Test] contact pt =',contacts(i,2,2:4)
        ! write(*,*)'[Bending Test] contact case =',contacts(i,1,9)

        
        !Case 2: Touching two fibers

        !Case 3: Touching three fibers
        
        
        end if
         
        !write(*,*)'Fiber bending_test',fiber_path(i,1,2)
       
       
    end subroutine bending_test


    subroutine fiber_mapping(i,dist,endpoint1,endpoint2,cog,p1,p2)!,bent)
        implicit none
        real(kind=8),dimension(4) :: p1,p2,endpoint1,endpoint2,cog
        integer,intent(in):: i
        !real(kind=8),dimension(max_fibers,1)::bent
    
        !local parameters
        real(kind=8),parameter              :: PI=4.0*atan(1.d0)
        integer(ip)::num_pt,m
        logical                               ::debug,aligned,neg
        real(kind=8)                :: d,mag_unit_fiber,length,length_original,lenght_tol
        real(kind=8),dimension(4,4) :: T,T_inv,Rx, Rx_inv, Ry, Ry_inv,Rz,Rz_neg
        real(kind=8),dimension(3) :: unit_fiber,unit_PQ,unit_cyl
        real(kind=8),dimension(4) :: axis,test,pt_p1_neg,pt_p1_pos,pt_p2_neg,pt_p2_pos,temp
        real(wp),dimension(5)::fx,y,dist
        real(wp),dimension(0:40) :: x_new,f_new_default
        
        !real(kind=8),dimension(max_fibers,8)::fibers
        aligned=.false.
        neg=.false.
        debug=.false.
        num_pt=5
    
    
    if(debug)write(*,*)'endpoint1 =',endpoint1
    if(debug)write(*,*)'p1 =',p1
    if(debug)write(*,*)'cog =',cog
    if(debug)write(*,*)'p2 =',p2
    if(debug)write(*,*)'endpoint2',endpoint2
         
       
    length_original=norm2((endpoint1(1:3)-endpoint2(1:3)))
    lenght_tol=0.0001
    if(debug)write(*,*)'length =',length_original
    !Define the translation and inverse translation 
    T=0.0
    T(1,1)=1.d0
    T(1,4)=-endpoint1(1)
    T(2,4)=-endpoint1(2)
    T(3,4)=-endpoint1(3)
    T(2,2)=1.d0
    T(3,3)=1.d0
    T(4,4)=1.d0
    
    T_inv=0.0
    T_inv(1,1)=1.d0
    T_inv(1,4)=endpoint1(1)
    T_inv(2,4)=endpoint1(2)
    T_inv(3,4)=endpoint1(3)
    T_inv(2,2)=1.d0
    T_inv(3,3)=1.d0
    T_inv(4,4)=1.d0
    
    !Translate end points
    endpoint1=matmul(T,endpoint1)
    endpoint2=matmul(T,endpoint2)
    p1=matmul(T,p1)
    p2=matmul(T,p2)
    cog=matmul(T,cog)
    debug=.true.
    if(debug)write(*,*)'Translated Endpoint 1 =',endpoint1
    if(debug)write(*,*)'Translated Endpoint 2 =',endpoint2
    test=0.d0
    test(3)=-0.1
    test(4)=1.d0
    !Defining a unit vector at the contact point pointing in the direction 
! of the top  fiber

unit_fiber(1)=endpoint2(1)-endpoint1(1)
unit_fiber(2)=endpoint2(2)-endpoint1(2)
unit_fiber(3)=endpoint2(3)-endpoint1(3)
mag_unit_fiber=sqrt(unit_fiber(1)**2+unit_fiber(2)**2+unit_fiber(3)**2)
unit_fiber=unit_fiber/mag_unit_fiber
if(debug)write(*,*)'unit_fiber =',unit_fiber

!Calculating a unit vector at the contact point pointing along the axis of rotation

axis(1)= unit_fiber(2)  
axis(2)= -unit_fiber(1) 
axis(3)= 0.d0
axis(4)=1.d0

axis(1:3)=(axis(1:3)/norm2(axis(1:3)))
!if (debug) if(debug)write(*,*)'Axis',axis
d=sqrt(axis(2)**2+axis(3)**2)

if(debug)write(*,*)'d =',d
if(debug)write(*,*)'axis =',axis
!Defining the rotation matrix

Rx=0.d0
Rx(1,1)=1.d0
Rx(2,2)=axis(3)/d
Rx(2,3)=-axis(2)/d
Rx(3,2)=axis(2)/d
Rx(3,3)=axis(3)/d
Rx(4,4)=1.d0

Rx_inv=0
Rx_inv(1,1)=1
Rx_inv(2,2)=axis(3)/d
Rx_inv(2,3)=axis(2)/d
Rx_inv(3,2)=-axis(2)/d
Rx_inv(3,3)=axis(3)/d
Rx_inv(4,4)=1

Ry=0
Ry(1,1)=d
Ry(1,3)=-axis(1)
Ry(2,2)=1
Ry(3,1)=axis(1)
Ry(3,3)=d
Ry(4,4)=1


Ry_inv=0
Ry_inv(1,1)=d
Ry_inv(1,3)=axis(1)
Ry_inv(2,2)=1
Ry_inv(3,1)=-axis(1)
Ry_inv(3,3)=d
Ry_inv(4,4)=1

if (d.gt.0.00001) then
    endpoint1=matmul(Rx,endpoint1)
    endpoint2=matmul(Rx,endpoint2)
    p1=matmul(Rx,p1)
    p2=matmul(Rx,p2)
    cog=matmul(Rx,cog)
    test=matmul(Rx,test)
else

end if

if(debug)write(*,*)'Rot. x Endpoint 2 =',endpoint2

endpoint1=matmul(Ry,endpoint1)
endpoint2=matmul(Ry,endpoint2)
p1=matmul(Ry,p1)
p2=matmul(Ry,p2)
cog=matmul(Ry,cog)
test=matmul(Ry,test)
if(debug)write(*,*)'Rot. y unit =',endpoint2
if(debug)write(*,*)'TEST =',test

! This is where the bending happens
pt_p1_neg=1.d0
pt_p1_pos=1.d0
pt_p2_neg=1.d0
pt_p2_pos=1.d0
pt_p1_neg(1:3)=p1(1:3)-((p1(1:3)-endpoint1(1:3))/(norm2((p1(1:3)-endpoint1(1:3)))))
pt_p1_pos(1:3)=p1(1:3)+((p1(1:3)-endpoint1(1:3))/(norm2((p1(1:3)-endpoint1(1:3)))))
pt_p2_neg(1:3)=p2(1:3)+((p2(1:3)-endpoint2(1:3))/(norm2((p2(1:3)-endpoint2(1:3)))))
pt_p2_pos(1:3)=p2(1:3)-((p2(1:3)-endpoint2(1:3))/(norm2((p2(1:3)-endpoint2(1:3)))))
!Step one determine what axis the fiber is aligned with

if (abs(endpoint2(1)).lt.0.0001)then
        if(debug)write(*,*)'Aligned with y axis'
        aligned=.true.

end if
endpoint1(1:3)=endpoint1(1:3)+test(1:3)*dist(1)!**2.d0
p1(1:3)=p1(1:3)+test(1:3)*dist(2)!**2.d0
cog(1:3)=cog(1:3)+test(1:3)*dist(3)!**2.d0
p2(1:3)=p2(1:3)+test(1:3)*dist(4)!**2.d0
endpoint2(1:3)=endpoint2(1:3)+test(1:3)*dist(5)!**2.d0
if (debug)write(*,*)'Endpoint 2 = ',endpoint2

if (aligned)then
fx(1)=endpoint1(2)
fx(2)=p1(2)
fx(3)=cog(2)
fx(4)=p2(2)
fx(5)=endpoint2(2)

y(1)=endpoint1(1)
y(2)=p1(1)
y(3)=cog(1)
y(4)=p2(1)
y(5)=endpoint2(1)
write(*,*)'Here 1438'
else
    if(endpoint2(1).lt.0.d0)then
    neg=.true.
        fx(1)=-endpoint1(1)
        fx(2)=-p1(1)
        fx(3)=-cog(1)
        fx(4)=-p2(1)
        fx(5)=-endpoint2(1)
    else
        fx(1)=endpoint1(1)
        fx(2)=p1(1)
        fx(3)=cog(1)
        fx(4)=p2(1)
        fx(5)=endpoint2(1)
    end if
        y(1)=endpoint1(2)
        y(2)=p1(2)
        y(3)=cog(2)
        y(4)=p2(2)
        y(5)=endpoint2(2)
end if


if(debug)write(*,*)'>>>>>> X = ',fx
if(debug)write(*,*)'y = ',y
        call sean_tests(fx,y,num_pt,x_new,f_new_default)
        if(neg)then
            x_new=-x_new
            fx=-fx
        end if
temp=1.d0
debug=.false.
if (aligned)then
    do m=0,40
        if(debug)write(*,*)'X=',x_new(m)
        if(debug)write(*,*)'fx =',f_new_default(m)
        temp(2)=x_new(m)
        temp(1)=f_new_default(m)
        temp(3)=0.d0
        temp=matmul(Ry_inv,temp)
        if(d.gt.0.00001)temp=matmul(Rx_inv,temp)
        temp=matmul(T_inv,temp)
        fiber_path(i,m+1,1:3)=temp(1:3)
        if(debug)write(*,*)'Fiber path = ',fiber_path(i,m+1,1:3)
    end do
else
        do m=0,40
            if(debug)write(*,*)'X=',x_new(m)
            if(debug)write(*,*)'fx =',f_new_default(m)
           temp(1)=x_new(m)
           temp(2)=f_new_default(m)
           temp(3)=0.d0
           temp=matmul(Ry_inv,temp)
           if(d.gt.0.00001)temp=matmul(Rx_inv,temp)
           temp=matmul(T_inv,temp)
           fiber_path(i,m+1,1:3)=temp(1:3)
           if(debug)write(*,*)'Fiber path = ',fiber_path(i,m+1,1:3)
        end do
end if

endpoint1=matmul(Ry_inv,endpoint1)
if(debug)write(*,*)' =',axis
endpoint2=matmul(Ry_inv,endpoint2)
p1=matmul(Ry_inv,p1)
p2=matmul(Ry_inv,p2)
cog=matmul(Ry_inv,cog)
if(debug)write(*,*)'Inverse y Endpoint 2 =',endpoint2
if (d.gt.0.00001) then
    endpoint1=matmul(Rx_inv,endpoint1)
    endpoint2=matmul(Rx_inv,endpoint2)
    p1=matmul(Rx_inv,p1)
    p2=matmul(Rx_inv,p2)
    cog=matmul(Rx_inv,cog)
    if(debug)write(*,*)'Inverse x Endpoint 2 =',endpoint2
    
else
end if
endpoint1=matmul(T_inv,endpoint1)
endpoint2=matmul(T_inv,endpoint2)
p1=matmul(T_inv,p1)
p2=matmul(T_inv,p2)
cog=matmul(T_inv,cog)
if(debug)write(*,*)'Inverse Translate Endpoint 2 =',endpoint2
length=norm2(endpoint2(1:3)-endpoint1(1:3))
if(debug)write(*,*)'length =',length
if(debug)write(*,*)'endpoint1 =',endpoint1
if(debug)write(*,*)'p1 =',p1
if(debug)write(*,*)'cog =',cog
if(debug)write(*,*)'p2 =',p2
if(debug)write(*,*)'endpoint2',endpoint2

fibers(i,1:3)=endpoint1(1:3)
fibers(i,5:7)=endpoint2(1:3)
control_points(i,1,1:3)=endpoint1(1:3)
control_points(i,2,1:3)=p1(1:3)
control_points(i,3,1:3)=cog(1:3)
control_points(i,4,1:3)=p2(1:3)
control_points(i,5,1:3)=endpoint2(1:3)
write(*,*)control_points(i,1,1:3)

!write(*,*)'Fiber >>>>',fiber_path(i,1,2)
! 210 format('_test ',f0.3,',',f0.3,',',f0.3,' ',f0.3,',',f0.3,',',f0.3,' '&
!            ,f0.3,',',f0.3,',',f0.3,' ',f0.3,',',f0.3,',',f0.3,' ',f0.3,','&
!            ,f0.3,',',f0.3,' ',f0.3,',',f0.3,',',f0.3,' ',f0.3,',',f0.3,','&
!            ,f0.3,' ',f0.3,',',f0.3,',',f0.3,' ',f0.3,',',f0.3,',',f0.3,' '&
!            ,f0.3,',',f0.3,',',f0.3,' ',f0.3,',',f0.3,',',f0.3,' ',f0.3)

    ! open(unit=16,file="test.scr", status='unknown')
    ! write(16,201)endpoint1(1:3),p1(1:3),cog(1:3),p2(1:3),endpoint2(1:3),fibers(i,4)
    ! close(16)
    
end subroutine fiber_mapping

subroutine  sean_tests(fx,y,num_pt,x_new,f_new_default)
implicit none

    integer(ip),intent(in)::num_pt
    real(wp),dimension(num_pt),intent(in)::fx,y

    integer :: i    !! counter

     integer(ip),parameter :: kx  = 3    !! x bspline order
     integer(ip),parameter :: nx= 5      !! number of points in x dimension
     real(wp),dimension(num_pt):: x       !! [0,20,40,60,80,100]
     !real(wp),dimension(num_pt)::fx,y
     real(wp),dimension(num_pt)    :: fcn
     type(bspline_1d)          :: s_default
     real(wp),dimension(0:40),intent(out) :: x_new,f_new_default
     !,f_actual
     real(wp)                  :: xval,step_size
     integer(ip)               :: iflag
     !write(*,*)'Here 1'
     x=fx
     fcn = y
     step_size=(x(num_pt)-x(1))/40
     !write(*,*)x
     !write(*,*)'x =',x
     !write(*,*)'y =',fcn
     !write(*,*)'Here 2'
     call s_default%initialize(x,fcn,kx,iflag)  !default (not-a-knot)
     if(iflag/=0)write(*,*)'iflag =', iflag
     if (iflag/=0) error stop 'error initializing s_default'
     !write(*,*)'Here 3'
     do i=0,40

        xval     = x(1)+i*step_size
        x_new(i) = xval

        !f_actual(i) = test_func(xval)

        call s_default%evaluate(xval,0_ip,f_new_default(i),iflag)
        if (iflag/=0) error stop 'error evaluating s_default'

     end do


end subroutine sean_tests

subroutine fiber_mapping_alt(i,dist,endpoint1,endpoint2,cog,p1,p2)!,bent)
    implicit none
    real(kind=8),dimension(4) :: p1,p2,endpoint1,endpoint2,cog
    integer,intent(in):: i
    !real(kind=8),dimension(max_fibers,1)::bent

    !local parameters
    real(kind=8),parameter              :: PI=4.0*atan(1.d0)
    integer(ip)::num_pt,m
    logical                               ::debug,aligned,neg
    real(kind=8)                :: d,mag_unit_fiber,length,length_original,lenght_tol
    real(kind=8)                ::alpha,theta
    real(kind=8),dimension(4,4) :: T,T_inv
    real(kind=8),dimension(4,4) ::ba, ba_trans
    real(kind=8),dimension(3) :: unit_fiber,unit_PQ,unit_cyl
    real(kind=8),dimension(4) :: axis,test,pt_p1_neg,pt_p1_pos,pt_p2_neg,pt_p2_pos,temp
    real(wp),dimension(5)::fx,y,dist
    real(wp),dimension(0:40) :: x_new,f_new_default

    num_pt=5
    test=0.d0
    test(3)=-0.1
    test(4)=1.d0
    debug=.true.
    neg=.false.
    
    if(debug)write(*,*)'endpoint1 =',endpoint1
    if(debug)write(*,*)'p1 =',p1
    if(debug)write(*,*)'cog =',cog
    if(debug)write(*,*)'p2 =',p2
    if(debug)write(*,*)'endpoint2',endpoint2
         
       
    length_original=norm2((endpoint1(1:3)-endpoint2(1:3)))
    lenght_tol=0.0001
    if(debug)write(*,*)'length =',length_original
    !Define the translation and inverse translation 
    T=0.0
    T(1,1)=1.d0
    T(1,4)=-endpoint1(1)
    T(2,4)=-endpoint1(2)
    T(3,4)=-endpoint1(3)
    T(2,2)=1.d0
    T(3,3)=1.d0
    T(4,4)=1.d0
    
    T_inv=0.0
    T_inv(1,1)=1.d0
    T_inv(1,4)=endpoint1(1)
    T_inv(2,4)=endpoint1(2)
    T_inv(3,4)=endpoint1(3)
    T_inv(2,2)=1.d0
    T_inv(3,3)=1.d0
    T_inv(4,4)=1.d0
    
    !Translate end points
    
    endpoint1=matmul(T,endpoint1)
    endpoint2=matmul(T,endpoint2)
    p1=matmul(T,p1)
    p2=matmul(T,p2)
    cog=matmul(T,cog)

    theta=atan(-endpoint2(3)/endpoint2(2))

    alpha=atan(-((cos(theta)*endpoint2(2)-sin(theta)*endpoint2(3))/endpoint2(1)))
    ba=0.d0
    ba(1,1)=cos(alpha)
    ba(1,2)=-sin(alpha)*cos(theta)
    ba(1,3)=sin(alpha)*sin(theta)
    ba(2,1)=sin(alpha)
    ba(2,2)=cos(alpha)*cos(theta)
    ba(2,3)=-cos(alpha)*sin(theta)
    
    ba(3,2)=sin(theta)
    ba(3,3)=cos(theta)
    ba(4,4)=1.d0
    ba_trans=transpose(ba)
    if (endpoint2(2).gt.0.0.and.endpoint2(2).gt.0.0)then
        endpoint1=matmul(ba,endpoint1)
        endpoint2=matmul(ba,endpoint2)
        p1=matmul(ba,p1)
        p2=matmul(ba,p2)
        cog=matmul(ba,cog)
        test=matmul(ba,test)
    end if
    if(endpoint2(1).lt.0)neg=.true.

    endpoint1(1:3)=endpoint1(1:3)+test(1:3)*dist(1)!**2.d0
    p1(1:3)=p1(1:3)+test(1:3)*dist(2)!**2.d0
    cog(1:3)=cog(1:3)+test(1:3)*dist(3)!**2.d0
    p2(1:3)=p2(1:3)+test(1:3)*dist(4)!**2.d0
    endpoint2(1:3)=endpoint2(1:3)+test(1:3)*dist(5)!**2.d0
    if(neg) then
        fx(1)=-endpoint1(1)
        fx(2)=-p1(1)
        fx(3)=-cog(1)
        fx(4)=-p2(1)
        fx(5)=-endpoint2(1)

        y(1)=endpoint1(2)
        y(2)=p1(2)
        y(3)=cog(2)
        y(4)=p2(2)
        y(5)=endpoint2(2)
    else
        fx(1)=endpoint1(1)
        fx(2)=p1(1)
        fx(3)=cog(1)
        fx(4)=p2(1)
        fx(5)=endpoint2(1)

        y(1)=endpoint1(2)
        y(2)=p1(2)
        y(3)=cog(2)
        y(4)=p2(2)
        y(5)=endpoint2(2)
    end if

    if(debug)write(*,*)'>>>>>> X = ',fx
    if(debug)write(*,*)'y = ',y
    call sean_tests(fx,y,num_pt,x_new,f_new_default)
    temp=1.d0
    do m=0,40
        if(debug)write(*,*)'X=',x_new(m)
        if(debug)write(*,*)'fx =',f_new_default(m)
        if (neg) then
            temp(1)=-x_new(m)
        else
            temp(1)=x_new(m)
        end if
        temp(2)=f_new_default(m)
        temp(3)=0.d0
        if (endpoint2(2).gt.0.0.and.endpoint2(2).gt.0.0)then
            temp=matmul(ba_trans,temp)
        end if
        temp=matmul(T_inv,temp)
        fiber_path(i,m+1,1:3)=temp(1:3)
        if(debug)write(*,*)'Fiber path = ',fiber_path(i,m+1,1:3)
    end do
    if (endpoint2(2).gt.0.0.and.endpoint2(2).gt.0.0)then
        endpoint1=matmul(ba_trans,endpoint1)
        endpoint2=matmul(ba_trans,endpoint2)
        p1=matmul(ba_trans,p1)
        p2=matmul(ba_trans,p2)
        cog=matmul(ba_trans,cog)
    end if
    
    endpoint1=matmul(T_inv,endpoint1)
    endpoint2=matmul(T_inv,endpoint2)
    p1=matmul(T_inv,p1)
    p2=matmul(T_inv,p2)
    cog=matmul(T_inv,cog)
    
    if(debug)write(*,*)'endpoint1 =',endpoint1
    if(debug)write(*,*)'P 1 =',p1
    if(debug)write(*,*)'P 2 =',p2
    if(debug)write(*,*)'COG =',cog
    if(debug)write(*,*)'endpoint2',endpoint2
    debug=.false.

    fibers(i,1:3)=endpoint1(1:3)
    fibers(i,5:7)=endpoint2(1:3)
    control_points(i,1,1:3)=endpoint1(1:3)
    control_points(i,2,1:3)=p1(1:3)
    control_points(i,3,1:3)=cog(1:3)
    control_points(i,4,1:3)=p2(1:3)
    control_points(i,5,1:3)=endpoint2(1:3)



end subroutine fiber_mapping_alt

end module fibers_place

