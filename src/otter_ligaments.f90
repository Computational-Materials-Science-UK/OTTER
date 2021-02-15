!!! OTTER - otter_ligaments.f90
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
! otter_spheres.f90 - module providing otter_ligaments structure generation version
!   v0.1 - Original OTTER variant, with thanks to many...
!
!   vb.1 - Stub for integrated OTTER
!   vb.2 - Beck completed integration... 02/15/2021

! External Dependencies:
!	None
! Internal Dependencies:
!	otter_main.f90

module ligaments_place

    use otter_globals
    use otter_input
    use otter_ligaments_globals

    contains

    !!!!!
    ! otter_ligaments - main subroutine for OTTER version dropping cylindrical ligaments
    !   v0.1 - Original OTTER variant, with thanks to many...
    !        1 : random spheres with buffer
    !        2 : added truncated cones...
    !        3 : added script output for AutoCAD.. 02/08/17
    !        5 : plate removal for displacement boundary conditions, 3/29/2017
    !        6 : changed # of nn's to store to depend on CN... 7/25/2017
    !	vb.1 - Working version - Beck re-work to integrate in main OTTER. 2/14/2021

    ! External Dependencies:
    !	none
    ! Internal Dependencies:
    !	[TBA]
    subroutine otter_ligaments()

        implicit none 

        !!!! Calling parameters
        !!!! Constants
        !!!! Variables
        real(kind=DBL)              :: sphere_range, &
                                       n_c_star
        integer                     :: nn_num, num_bad, &
                                       rve, lig_count, end_sphere, curr_sphere, &
                                       n_c_tot, n_sphere_in, &
                                       i,j,k,l,jj, &
                                       scr_unit,dat_unit
        character(len=200)          :: full_scr_name,full_nnc_name,saveas,exportstl,full_name
        logical                     :: good, in_box, found
        character                   :: cont, num_char
        character(len=50)           ::  date_time

        real(kind=8),allocatable,dimension(:,:)     :: tempnn

        real(kind=8),dimension(4)                   :: dist
        real(kind=8),dimension(2,3)                 :: bigbox
        real(kind=8),dimension(3)                   :: shift
        

        !!!! Begin code... !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        write (out_unit,*) ' Welcome to Ligaments vb.1!'
        write (out_unit,*) ' '

        !!!! Input !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Get Global/General input, see otter_input
        call get_gen_input()
        ! Get spheres input..., see otter_input
        call get_ligaments_input()

        !!!! Prep tasks... !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! normalize sphere numbers to preserve sphere density in big box
        sphere_num=NINT(box_vol_factor*sphere_num)
        sphere_range=max_rad-min_rad
        shift(:)=box_length(:)*((box_vol_factor**(1.d0/3.d0))-1.d0)
        box_length(:)=box_length(:)+shift(:)
        shift(:)=-0.5d0*shift(:)

        nn_num=2*num_ligs

        ! Output inputs...


        !!!! Allocate Arrays... !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! spheres:  sphere data for each of sphere_num spheres: x-pos, y-pos, z-pos, radius
        allocate(spheres(sphere_num,4))

        ! nn:  nearest neighbor table contain nn_num neighbors for each of sphere_num spheres:
        !      nn atom number, x-vec, y-vec, z-vec, distance
        allocate(nn(sphere_num,nn_num,5))

        ! tempnn:  temp arrary for nn table data while determining if sphere is good:
        !          nn atom number, x-vec, y-vec, z-vec, distance
        allocate(tempnn(sphere_num,5))

        ! ligmap:  table of which of a sphere's nn_num near neighbors are connected to the sphere
        !  For each of sphere_num spheres the first nn_num fields are for each of the first nn_num neighbors. The fields are 0 if not yet connected, 1 if connected. The nn_num+1st field is the total number of ligaments connected to the sphere, and may be higher than the sum of the other fields of other spheres are connected (should be rare).
        allocate(lig_map(sphere_num,nn_num+1))

        ! ligaments:  analog to spheres... details describing each ligament, of where there are a maximum of neigh_model_num*sphere_num
        !       each ligament has 2 end points, each end point has x-pos, y-pos, z-pos, radius
        allocate(ligaments(num_ligs*sphere_num,2,4))

        !!!! Generate RVEs
        ! Generate structures...
        write (out_unit,*) ' Generating Structure(s)...'

        ! Set-up output for Auto-CAD multi-script
        ! **** THERE is a weakness here... output should be centralized, and it should check for existing files before crashing upon discovery....
        scr_unit=53
        dat_unit=54
        full_scr_name=trim(rves_batch_name)//".scr"
        full_nnc_name=trim(rves_batch_name)//".out"
        open(unit=scr_unit,file=trim(full_scr_name),status='new',action='write')
        open(unit=dat_unit,file=trim(full_nnc_name),status='new',action='write')

        ! Write dat file header...
        call fdate(date_time)
        write(dat_unit,*) rves_batch_name
        write(dat_unit,*) ' Date generated: ',date_time
        write(dat_unit,*) ' Number of RVEs: ',rves_num
        write(dat_unit,*) ' Box length(x,y,z), Big Box Volume factor:'
        write(dat_unit,*) box_length(:),box_vol_factor
        write(dat_unit,*) ' Sphere num, sphere buffer, sphere inflation:'
        write(dat_unit,*) sphere_num, sphere_buffer, sphere_inflation
        write(dat_unit,*) ' Min rad, max rad:'
        write(dat_unit,*) min_rad, max_rad
        write(dat_unit,*) ''
        write(dat_unit,*) ' RVE data: RVE #, n_c_tot, n_sphere_in, n_c_star '

        !Begin building structures, one model at a time...
        do rve=1,rves_num
            write(out_unit,*) 'RVE #',rve,'...'
            
            ! Set-up output for individual RVE
            write(num_char,"(i4)") rve ! cast a four character version of the rve number
            full_name=trim(rves_batch_name)//trim(adjustl(num_char))

            ! Initialize for each structure
            nn(:,:,:)=0.0
            nn(:,:,5)=999999999.9  ! Prep columns used for sorting distances
            i=1
            num_bad=0

            !!!! Construct Spheres !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do while (i .le. sphere_num)  ! Create spheres until sphere_num are created...
                tempnn(:,:)=0.0
                tempnn(:,5)=9999999999.9 ! fails to sort properly if box_length is VERY large

                spheres(i,1)=rand()*box_length(1)+shift(1)
                spheres(i,2)=rand()*box_length(2)+shift(2)
                spheres(i,3)=rand()*box_length(3)+shift(3)
                spheres(i,4)=min_rad+rand()*sphere_range
                if(debug) write(*,*) 'Sphere loop: ',i,spheres(i,:)

                ! Determine if potential sphere i overlaps with any sphere j =>bad
                ! ...also store all nearest neighbors of sphere i in tempnn!...
                good=.true.   !assume new sphere is good
                j=1
                tempnn(:,5)=999999999.9

                do while (good .and. (j .lt. i))  !only need to compare to previous spheres
                    dist(1:3)=spheres(i,1:3)-spheres(j,1:3)
                    dist(4)=sqrt(dist(1)**2+dist(2)**2+dist(3)**2)-spheres(i,4)-spheres(j,4)
                    if(debug) write(*,*) 'NN loop: ',j,dist(:)
                    !Check if sphere being considered (i) is too close to test sphere (j)
                    if (dist(4) .lt. sphere_buffer) then  
                        good=.false.
                        if(debug) write(*,*) '    NN -> not good sphere'
                        num_bad=num_bad+1
                    else  ! If i is not too close to j, add j to the tempnn list in it's position...
                        k=i   ! the farthest possible neighber j can be is the ith, as we are adding the ith sphere. Now we'll count down looking for when the jth sphere is no longer closer than the k-1st neighbor, so the jth sphere is the ith sphere's kth nearest neighbor 
                        do while ((k .gt. 1) .and. (tempnn(k-1,5) .gt. dist(4)))
                            tempnn(k,:)=tempnn(k-1,:)  ! if it's not the kth nn, then the k-1st must become the  kth
                            k=k-1
                        end do
                        tempnn(k,1)=j ! now we're out of the while, j is the k nn of i
                        tempnn(k,2:5)=dist(:)
                        if(debug) write(*,*) '    NN -> ',j,' is NN',k,':',tempnn(k,:)
                        if(debug) read(*,*)
                    endif
                    j=j+1  ! counting up through spheres previously added checking if i overlaps j
                end do

                ! If it's still good, we need to add the sphere and it's nn's to the lists...
                ! The hard part is that we have to reciprically update the already existing spheres nn lists.
                if (good) then
                    nn(i,1:nn_num,:)=tempnn(1:nn_num,:)  ! store the new ith sphere's nn_num  nn's in the table.... 
                    !if (debug) write(*,'(i3,a4)',ADVANCE='NO') i,'... '
                    if (debug) write(*,*) '    Adding good sphere...',nn(i,1,:)

                    do j=1,i-1
                        ! starting with the ith spheres 1st nn, we'll get the sphere's id. Then, counting down as we did above, we'll insert the new ith sphere in it's place in the full nn table for the "curr_sphere"... note that we start counting down with k-1=nn_num, because we only are keeping nn_num spheres in the neighbor list.
                        k=nn_num+1
                        curr_sphere=tempnn(j,1) 
                        do while ((k .gt. 1) .and. (nn(curr_sphere,k-1,5) .gt. tempnn(j,5)))
                            if (k .lt. nn_num+1) nn(curr_sphere,k,:)=nn(curr_sphere,k-1,:)
                            k=k-1
                        enddo
                        if (k .lt. nn_num+1) then ! w/ above, only shift+insert if i is in curr_sphere's nn_num nn's
                            nn(curr_sphere,k,1)=i
                            nn(curr_sphere,k,2:4)=0-tempnn(j,2:4)
                            nn(curr_sphere,k,5)=tempnn(j,5)
                        endif
                    end do

                    i=i+1  ! increment the number of good spheres
                    num_bad=0  ! reset the num_bad count
                end if  ! ... end of the "if good" loop

                if (num_bad .ge. 100000) then
                    write(*,*) ' 100000 consecutive bad spheres!... Continue?'
                    read(*,*) cont
                    if ((cont .eq. 'y') .or. (cont .eq. 'Y')) then
                        num_bad=0
                    else
                        exit  ! exits the create sphere loop if user says so...
                    end if
                end if ! end if too many bads in a row...
            end do  ! End sphere i generation loop for structure rve

            !!!! Construct Cones !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            write(*,*)
            if (num_bad .ne. 100000) then ! only generate cones if user hasn't exitted....
      
                write(*,*) ' Generating cones...'
                lig_map(:,:)=0
                lig_count=0
                ligaments(:,:,:)=-1.0
                do i=1,num_ligs
                    ! cycle spheres, do 1 lig for every sphere, then another, up to num_ligs...
                    do j=1,sphere_num
                        ! Check if sphere is in final box
                        in_box=.true.
                        do jj=1,3
                            if ((spheres(j,jj) .lt. 0.0) .or. (spheres(j,jj) .gt. box_length(j)+2.d0*shift(j))) in_box=.false.
                        end do
                        ! does the sphere already have the number of ligaments for this round?
                        if ((lig_map(j,nn_num+1) .lt. i) .and. in_box) then
                            ! Find the first unconnected NN... max is nn_num for how many nn's...
                            k=1
                            do while ((k .lt. nn_num+1) .and. (lig_map(j,k) .ne. 0)) ! searching for 1st unconnected k
                                k=k+1
                            end do
                            ! if there is an unconnected ligament, add it...
                            if (k .lt. nn_num+1) then
                                lig_count=lig_count+1
                                ligaments(lig_count,1,1:3)=spheres(j,1:3)
                                ligaments(lig_count,1,4)=spheres(j,4)/(1+sphere_inflation)
                                end_sphere=nn(j,k,1)
                                ligaments(lig_count,2,1:3)=spheres(end_sphere,1:3)
                                ligaments(lig_count,2,4)=spheres(end_sphere,4)/(1+sphere_inflation)
                                lig_map(j,k)=1
                                lig_map(j,nn_num+1)=lig_map(j,nn_num+1)+1  !increment total ligaments on node
                                ! now mark off the other end of the ligament in the lig_map
                                lig_map(end_sphere,nn_num+1)=lig_map(end_sphere,nn_num+1)+1
                                found=.false.
                                l=1
                                do while ((l .lt. nn_num+1) .and. (.not. found))
                                    if (nn(end_sphere,l,1) .eq. j) then
                                        found=.true.
                                        lig_map(end_sphere,l)=1
                                    else
                                        l=l+1
                                    endif
                                enddo ! find other end of lig...
                            endif ! check if available unconnected lig
                        endif ! stuff to do if sphere doesn't have i ligs and is in box
                    enddo ! add lig i to each sphere j
                enddo ! construct i out of num_lig ligs for each node...

            end if ! construct ligs for structure rve unless user says to stop.

            !!!! Write Output !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

202         format('_sphere ',f0.3,',',f0.3,',',f0.3,' ',f0.3)
203         format('_cone ',f0.3,',',f0.3,',',f0.3,' ',f0.3,' T ',f0.3,' A ',f0.3,',',f0.3,',',f0.3)
!204         format('_box ',f0.3,',',f0.3,',',f0.3,' ',f0.3,',',f0.3,',',f0.3)
            write(scr_unit,'(a)') '-osnap off'
            write(scr_unit,'(a)') '_zoom -300,-300,-300 500,500,500' 
            write(scr_unit,'(a)') '_FACETRES 6'
            do i=1,sphere_num  
                in_box=.true.
                do jj=1,3
                    if ((spheres(i,jj) .lt. (0.0-spheres(i,4))) .or. (spheres(i,jj) .gt. &
                        (box_length(j)+2.d0*shift(j)+spheres(i,4)))) in_box=.false.
                end do
                if (in_box) write(scr_unit,202) spheres(i,1:4)
            end do
            write(scr_unit,'(a)') '_union all '
            do j=1,lig_count
                write(scr_unit,203) ligaments(j,1,1:4),ligaments(j,2,4),ligaments(j,2,1:3)
                if (mod(j,20) .eq. 0) then
                    write(scr_unit,'(a)') '_union all '
                end if
            end do

            write(scr_unit,'(a)') '_union ALL '
            write(scr_unit,'(a)') 'FILEDIA 0'
            saveas='_saveas  '//trim(out_path)//trim(adjustl(full_name))//'.dwg'
            ! write(out_unit,'(a)') trim(adjustl(saveas))

            write(scr_unit,'(a)') '_group create network  all '
            do i=1,3
                bigbox(1,i)=0.d0+shift(i)
                bigbox(2,i)=box_length(i)+shift(i)
            end do

301         format('_box ',f0.4,',',f0.4,',',f0.4,' ',f0.4,',',f0.4,',',f0.4)
            write(scr_unit,301) bigbox(1,1),bigbox(1,2),bigbox(1,3),bigbox(2,1),bigbox(2,2),bigbox(2,3)
            write(scr_unit,'(a)') '_group create bigbox  last '
            write(scr_unit,301) 0.,0.,0.,box_length(1)+2.d0*shift(1),box_length(2)+2.d0*shift(2),box_length(3)+2.d0*shift(3)
            write(scr_unit,'(a)') '_group create smallbox  last '
            write(scr_unit,'(a)') '_subtract g bigbox  g smallbox '
            write(scr_unit,'(a)') '_subtract g network  g bigbox '
            saveas='_saveas  '//trim(out_path)//trim(adjustl(full_name))//'np.dwg'
            exportstl='_export '//trim(out_path)//trim(adjustl(full_name))//'np.stl all'
            ! write(scr_unit,'(a)') '_move all  0,0,0 500,500,500'
            ! write(scr_unit,'(a)') trim(saveas)
            write(scr_unit,'(a,a)') trim(exportstl),' '

            write(scr_unit,'(a)') '_erase all '
            write(scr_unit,'(a)') '-purge g network y y'
            write(scr_unit,'(a)') '-purge g bigbox y y'
            write(scr_unit,'(a)') '-purge g smallbox y y'

            !!!! Collect RVE statistics

            n_c_tot=0
            n_sphere_in=0
            do i=1,sphere_num
                ! Check if sphere was in box for ligament creation
                in_box=.true.
                do jj=1,3
                    if ((spheres(j,jj) .lt. 0.0) .or. (spheres(j,jj) .gt. box_length(j)+2.d0*shift(j))) in_box=.false.
                end do
                if (in_box) then
                    n_c_tot=n_c_tot+int(lig_map(i,nn_num+1))
                    n_sphere_in=n_sphere_in+1
                end if
            end do
            n_c_star=n_c_tot/n_sphere_in

            write(dat_unit,*) rve, n_c_tot, n_sphere_in, n_c_star

        end do ! done with rve out of num_rves

        write(scr_unit,'(a)') 'FILEDIA 1'
        write(scr_unit,'(a)') '_close'

        
    end subroutine otter_ligaments

end module ligaments_place