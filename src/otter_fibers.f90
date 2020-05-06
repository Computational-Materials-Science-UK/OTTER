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
! otter_spheres.f90 - module providing otter_fibers structure generation version
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

        !!!! Begin code...

        write (out_unit,*) ' Welcome to Fibers vb.1!'
        write (out_unit,*) ' '

        ! Get Global/General input, see otter_input
        call get_gen_input()

        ! Get spheres input...
        write(out_unit,*) ' Minimum and maximum fiber radii? '
        read(in_unit,*) min_rad, max_rad
        write(out_unit,*) ' Minimum and maximum fiber length? '
        read(in_unit,*) min_length, max_length
        write(out_unit,*) ' Minimum overlap distance? '
        read(in_unit,*) min_olp
        write(out_unit,*) ' Step size for moving fibers? '
        read(in_unit,*) step_size
        write(out_unit,*) ' Boundary thicknesses (+x, -x, +y, -y, +z, -z): '
        read(in_unit,*) box_boundary(1:6)

        ! compute 'big box' lengths, volume; shift origin of cut box to (0,0,0)
        big_box_vol=1.d0
        do i=1,3
            box_range(i)=box_length+box_boundary(i*2-1)+box_boundary(i*2)
            big_box_vol=big_box_vol*box_range(i)
            shift(i)=0.d0-box_boundary(i*2)
        end do
        ! Estimate max number of fibers as big_box_vol/cylindrical plate containing smallest cylinder... no idea if this is any good...
        max_fibers=3*int(big_box_vol/( PI*(min_length/2.)**2.*(2.*min_rad) ))
        write(out_unit,*) ''
        write(out_unit,*) ' Big box lengths, volume: ',box_range(:),big_box_vol
        if (debug) write(out_unit,*) ' (debug) shift, max_fibers: ',shift(:),max_fibers

        ! Generate structures...
        write (out_unit,*) ' Generating Structure(s)...'

        !!! Allocate arrays:
        ! spheres:  sphere data for each of max_fibers spheres - nn_count is number of contacting spheres...
        !       x-pos, y-pos, z-pos, radius, nn_count, in_box (1=yes, 0=no)
        allocate(fibers(max_fibers,6))

        ! nn:  nearest neighbor table for each of max_fibers shperes:
        !   contact sphere number, distance, overlap distance
        ! Each sphere in theory has 3 nn bc it falls until it touches 3..., 
        ! NN are only detailed here for the sphere landing, though all contacts are counted in spheres(i,5)
        allocate(nn(max_fibers,5,3))

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
            fibers(:,6)=1.d0
            nn(:,:,:)=0.d0
            i=1
            num_bad=0

            radius_range=max_rad-min_rad
            length_range=max_length-min_length
            if(debug) write(*,*) 'otter_fibers: fiber_range,length_range: ',fiber_range,length_range

            ! Main loop... create fibers until fiber maximum is exceeded or until the last 0.05*max_fibers have been bad, that is, above the big_box.
            do while ((i .le. max_fibers) .and. (num_bad .le. int(0.05*max_fibers)))
                if (debug) write(out_unit,'(i3,a4)',ADVANCE='NO') i,'... '

                ! randomly place a fiber above box
                fibers(i,1)=rand()*box_range(1)+shift(1)
                fibers(i,2)=rand()*box_range(2)+shift(2)
                fibers(i,3)=box_range(3)+2.*max_rad+shift(3)
                fibers(i,4)=min_rad+rand()*radius_range
                fibers(i,9)=min_length+rand()*length_range
                dir_rand=rand()*PI
                fibers(i,5)=fibers(i,5)+fibers(i,9)*cos(dir_rand)
                fibers(i,6)=fibers(i,6)+fibers(i,9)*sin(dir_rand)
                fibers(i,7)=fibers(i,3)
                fibers(i,8)=min_rad+rand()*radius_range

                ! Drop fiber until contact, or hits bottom.
                do while (.not. stopped)

                    ! Check for contacts
                    call fiber_contact(i,min_olp)
                    if (debug) write(*,*) 'Re-entered program spheres_place from fiber_contact'

                    select case (num_nn(i))
                    case (0)
                        ! Drop fiber by step_size
                        fibers(i,3)=fibers(i,3)-step_size
                        fibers(i,7)=fibers(i,7)-step_size
                    case (1)
                        fibers(i,:)=rotate_fiber(i,step_size)
                    case (2)
                        fibers(i,:)=rotate_fiber(i,step_size)
                    case default
                        ! Kludge here... in theory, there could be three contacts and fiber could still roll... this should be rare...
                        if (debug) write(out_unit,*) ' Stopped...'
                        stopped=.true.
                    end select



                end do ! dropping loop

            end do ! fiber generation loop

        end do ! rves loop
                

    end subroutine otter_fibers

    function fiber_rotate(i,step_size) return(new_pos)
        implicit None

        use otter_globals
        use otter_fibers_globals
        use otter_math




    end function fiber_rotate


    subroutine fiber_contact(i,min_olp,contact,contact_pt3D)
        implicit none 

        use otter_globals
        use otter_fibers_globals
        use otter_math

        num_nn(i)=0
        nn_list(i,:)=0

        ! Check gross contact
        do j=1,i-1
            contact=.true.

            do k=1,3
                if ((fibers(i,k+4) .lt. fibers(j,k)-2.*max_rad) .or. (fibers(i,k) .gt. fibers(j,k+4)+2.*max_rad)) contact=.false.
            end do

            ! check fine contact...
            if (contact) then
                contact=.false.
                ! rounded end approximation...
                call segments_dist_3d(fibers(i,1:3),fibers(i,5:7),fibers(j,1:3),fibers(j,5:7),dist,near_ptP(1), near_ptP(2))
                near_pt3D(1,:)=fibers(i,1:3)+near_ptP(1)*(fibers(i,5:7)-fibers(i,1:3))
                near_pt3D(2,:)=fibers(j,1:3)+near_ptP(2)*(fibers(j,5:7)-fibers(j,1:3))
                contact_vec(:)=near_pt3D(2,:)-near_pt3D(1,:)
                fiber_vec(1,:)=near_pt3D(1,:)-(fibers(i,5:7)+near_ptP(1)*(fibers(i,1:3)-fibers(i,5:7)))
                fiber_vec(2,:)=near_pt3D(2,:)-(fibers(j,5:7)+near_ptP(2)*(fibers(j,1:3)-fibers(j,5:7)))
                radius(1)=fibers(i,4)+near_ptP(1)*(fibers(i,8)-fibers(i,4))
                radius(2)=fibers(j,4)+near_ptP(2)*(fibers(j,8)-fibers(j,4))
                if ((near_ptP(1) .eq. 0.d0) .or. (near_ptP(1) .eq. 1.d0)) then
                    face_vec(1,:)=contact_vec(:)-dotn(contact_vec(:),fiber_vec(1,:)/norm2(fiber_vec(1,:)))
                    if ((near_ptP(2) .eq. 0.d0) .or. (near_ptP(2) .eq. 1.d0)) then
                        ! end-to-end : Check if proj of contact in each face is less than radii...
                        face_vec(2,:)=([0.d0,0.d0,0.d0]-contact_vec(:))-dotn(([0.d0,0.d0,0.d0]-contact_vec(:)),fiber_vec(2,:)/norm2(fiber_vec(2,:)))
                        ! kludge here... overlap is more than min_olp becuase it's applied to both face_vecs...
                        if ((norm2(face_vec(1,:)) .lt. radius(1)-min_olp) .and. (norm2(face_vec(2,:)) .lt. radius(2)-min_olp)) then
                            ! end-to-end contact
                            contact=.true.
                            contact_pt3D(:)=0.5*contact_vec(:)-dotn(0.5*contact_vec(:),fiber_vec(1,:)/norm2(fiber_vec(1,:)))+near_pt3D(1,:)
                        end if
                    else
                        ! end-to-side : Check if edge point of i is within radius of j.
                        ! There is a small kludge here... the contact point is set to the corner of the contacting cylinder...
                        contact_pt3D(:)=near_pt3D(1,:)+radius(1)*face_vec(1,:)/norm2(face_vec(1,:))
                        if (pt2line(fibers(j,1:3),fibers(j,5:7),contact_pt3D(:)) .lt. radius(2)-min_olp) then
                            ! end-to-side contact
                            contact=.true.
                        endif
                    end if

                else if ((contact_pt(2) .eq. 0.d0) .or. (contact_pt(2) .eq. 1.d0)) then
                    ! end-to-side : Check if edge point of j is within radius of i.
                        ! There is a small kludge here... the contact point is set to the corner of the contacted cylinder...
                    contact_pt3D(:)=near_pt3D(2,:)+radius(2)*face_vec(2,:)/norm2(face_vec(2,:))
                    if (pt2line(fibers(i,1:3),fibers(i,5:7),contact_pt3D(:)) .lt. radius(1)-min_olp) then
                        ! end-to-side contact
                        contact=.true.
                    endif
                else
                    ! side-to-side contact
                    if (dist .lt. radius(1)+radius(2)-min_olp) then
                        ! side-to-side contact
                        contact=.true.
                        contact_pt3D(:)=near_pt3D(1,:)+contact_vec(:)/dist*(radius(1)-0.5*(radius(1)+radius(2)-dist))
                    end if
                end if

            end if ! end fine check

            if (contact) then
                num_nn(i)=num_nn(i)+1
                nn_list(i,num_nn(i))=j
            end if

        end do

    end subroutine




end module fibers_place

