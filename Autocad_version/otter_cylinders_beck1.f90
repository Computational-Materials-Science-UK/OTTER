program cylinders_place

    use otter_input
    use otter_cyl_rotate
    use otter_cyl_contact
    implicit none

    !constants

    !variables
    character(len=10) :: geocode
    real(kind=8),allocatable,dimension(:,:)::cyls,rotation_count
    integer,allocatable,dimension(:,:)::nn_count
    real(kind=8),allocatable,dimension(:,:,:)::nn
    real(kind=8),allocatable,dimension(:,:)::spheres
    real(kind=8),dimension(2) ::bigbox,smallbox
    real(kind=8),dimension(3)::r,distance_cyls,unit_cyls
    real(kind=8),dimension(3,3)::pq
    logical::debug,in_box,contact,test,try,input,small_box,rotation_error
    integer::num_models,cylinder_num,out_unit=51,model,i,j,jj,indx(2),n,onenn,l,x,cont_i,m
    integer::periodic_i,periodic_j,periodic_added,preserve_cylinder_num
    real(kind=8)::box, box_input,min_rad,max_rad,buffer,min_length,max_length,radius_range,length_range
    real(kind=8) :: box_range,shift,cyldist,arclength,sphere_range,normalization,drop,average_rotation
    real(kind=8)::d,k,rot_test,slopeone,slopetwo,slopemax,z_rand,dir_rand,q_plane_dot,sum_rotation
    real(kind=8)::box_center,box_top_corner,box_bottom_corner,pi,length_cyls
    real(kind=8),dimension(3,2)::p,q
    real(kind=8),dimension(3,1)::pr,qr
    character(len=18)::out_name
    character(len=200)::path,full_out_name,full_scr_name,full_nnc_name,full_nncd_name,exportstl,full_name,saveas
    character (len=8) :: srf_char,num_char
    integer,dimension(8)::values
    character(len=1)::lig_edge

    !new_contact
    integer :: num_contact
    real(kind=8),dimension(2,3) :: contact_pt, curr_pt, i_pt
    integer,dimension(2) :: i_num

    !version history
    integer,parameter :: version=2
     pi=acos(-1.0)
    !1: randomly places cylinders outside a box and allows them to fall about eachother
    debug=.false.
    test=.true.
    try=.true.
    input=.false.
    small_box=.false.

    geocode='fibers-1'

    write(*,*) ' Welcome to Otter!'
    write(*,*) '   Starting up... '
    write(*,*) '   Geocode: ',geocode
    
    read(*,*)
    !Begin User input
   
    if (.not. test) then
        call get_input(geocode,num_models,out_name,path,box,min_rad,max_rad,buffer,min_length,max_length,cylinder_num,lig_edge)
    else
        
        num_models=100
        out_name='cut_box_28_'
        path='./'
        box=30
        min_rad=5.d-1
        max_rad=5.d-1
        buffer=min_rad/4.0
        min_length=10
        max_length=10
        cylinder_num=1800
        lig_edge='y'
    end if
    cylinder_num=cylinder_num*9 !periodic boundaries for edges rising as fast as center... see below.
    
    !ligamented edges:
    shift=0.0
    
    ! Ligamented edges have ligaments that are truncted at the faces.  This is done by
    ! making a box with 2x the volume, and therefore cube-root-of-two-times the side length.
    ! To keep the density of spheres the same, the number of spheres is doubled.
    ! The "shift" places the corner of the unit box at the origin.
    if ((lig_edge .eq. 'y') .or. (lig_edge .eq. 'Y')) then
        cylinder_num=2*cylinder_num
        box_range=1.26*box
        shift=0-0.13*box
    else
        box_range=box
        shift=0.0
    end if

    !! structure:
    
    write (*,*) ' Generating Structure(s)...'
    write (*,*) cylinder_num
    
    
    !allocate arrays
    !cyls: cylinder data for each of the cylinder_num cylinders
    !       data is stored for both circular endpoint centers
    !       a,b,c, raidus1, radius2,x, y, and z [coordinate points], magnitude
    allocate(cyls(cylinder_num,9))

    !spheres:
    allocate(spheres(cylinder_num,4))

    !  nn_count:
    allocate(nn_count(cylinder_num,2))
    nn_count=0
    allocate(rotation_count(cylinder_num,1))
    rotation_count=0
    !nn:  nearest neighbor table for each of sphere_num shperes
    !       nn_num neighbors are stored with
    !               atom_number, rx,ry,rz,PQx,PQy,PQz,cyls two z
    !each atom gets 2 nn bc it falls until it touches 2, then all nn are only accounted for once
    allocate(nn(cylinder_num,2,8))
   
    call date_and_time(values=values)
     !values=(/2020, 5, 4, -240, 13, 10, 35, 71/)
     !write(*,*)values
     !read(*,*)
    call srand(values(8))

    radius_range=max_rad-min_rad
    length_range=max_length-min_length
    if(debug) write(*,*) 'ranges: ',radius_range,length_range
    !Preserve original values 
    preserve_cylinder_num = cylinder_num
    box_input = box
    !Begin building structures, one model at a timee...
    do model=23,num_models
        write(*,*) 'Model #',model,'...'
        write(num_char,"(i4)") model !cast a four character version of the model number
        write(*,*)'num_char',num_char
        if (debug) write(*,*)'Number of cylinders', cylinder_num
        rotation_error=.false.
        !Full script (i.e generates one script with info for every model)
        
        full_name=trim(out_name)//trim(adjustl(num_char))
        full_scr_name=trim(full_name)//".scr"
        full_nnc_name=trim(full_name)//".txt"
        full_out_name=trim(full_name)//".out"
        full_nncd_name=trim(full_name)//".plain"
        open(unit=out_unit,file=trim(full_scr_name),status='new',action='write')
        open(unit=47,file=trim(full_nnc_name),status='new',action='write')
        open(unit=50,file=trim(full_out_name),status='new',action='write')
        open(unit=46,file=trim(full_nncd_name),status='new',action='write')
        if (input) open(unit=16,file='cylinders.out', status='old', action='read')
        


    !initialize arrays...
        nn_count(:,:)=0.0
        nn(:,:,:)=0.0
        cylinder_num = preserve_cylinder_num
        box=box_input
        i=1
        !cylinder_num=2000
        !cylinder_num=cylinder_num*9
        contact=.false.
        write(*,*)'Number o cylinders', cylinder_num
    !create cylinders until they've all been created
        do while(i .le. cylinder_num)
            write(*,'(i3,a4)',ADVANCE='NO') i,'... '
            !cyls(i,9)=max_length+10
            !slopeone=10000
            !slopetwo=1
            !z_rand=0.1*rand()-0.1
            !do while (cyls(i,9) .gt. max_length)
                !randomly place two end points within the box range, on top of the box
            if (input) then
110                 format (i4, 9(ES23.16,1x))
                read (16,110)m, cyls(i,1:8)
            else
                 if (test) then
!~                      if(i.eq.1) then
                       
!~                          cyls(i,1)=1.0
!~                          cyls(i,2)=10.0
!~                          cyls(i,3)=box_range+1+max_length
!~                          cyls(i,4)=min_rad+(rand()*radius_range)
!~                          cyls(i,5)=min_rad+(rand()*radius_range)
!~                          cyls(i,9)=rand()*length_range+min_length
!~                          dir_rand=rand()*2.d0*4.d0*atan(1.d0)
!~                          cyls(i,6)=cyls(i,1)+cyls(1,9)
!~                          cyls(i,7)=cyls(i,2)
!~                          cyls(i,8)=box_range+1+max_length
!~                      elseif(i.eq.10)then
!~                          cyls(i,1)=1.25
!~                          cyls(i,2)=9.75
!~                          cyls(i,3)=box_range+1+max_length
!~                          cyls(i,4)=min_rad+(rand()*radius_range)
!~                          cyls(i,5)=min_rad+(rand()*radius_range)
!~                          cyls(i,9)=rand()*length_range+min_length
!~                          dir_rand=rand()*2.d0*4.d0*atan(1.d0)
!~                          cyls(i,6)=11.0!45.3554
!~                          cyls(i,7)=9.85 !36.3554
!~                          cyls(i,8)=box_range+1+max_length
!~                  else
!~                      cyls(i,1)=0.0
!~                      cyls(i,2)=2.0
!~                      cyls(i,3)=box_range+1+max_length
!~                      cyls(i,4)=min_rad+(rand()*radius_range)
!~                      cyls(i,5)=min_rad+(rand()*radius_range)
!~                      cyls(i,9)=rand()*length_range+min_length
!~                      dir_rand=rand()*2.d0*4.d0*atan(1.d0)
!~                      cyls(i,6)=10.0
!~                      cyls(i,7)=2.0
!~                      cyls(i,8)=box_range+1+max_length
!~                  endif
!~                      write(50,*)i,cyls(i,1:8)
!~                  else
                    cyls(i,1)=rand()*box_range+shift
                    cyls(i,2)=rand()*box_range+shift
                    cyls(i,3)=box_range+1+(rand()-0.5)*max_length
                    cyls(i,4)=min_rad+(rand()*radius_range)
                    cyls(i,5)=min_rad+(rand()*radius_range)
                    cyls(i,9)=rand()*length_range+min_length
                    dir_rand=rand()*2.d0*4.d0*atan(1.d0)
                    cyls(i,6)=cyls(i,1)+cyls(i,9)*cos(dir_rand)
                    cyls(i,7)=cyls(i,2)+cyls(i,9)*sin(dir_rand)
                    cyls(i,8)=cyls(i,3)+cyls(i,9)*sin((rand()-0.5)*(12*pi/180))
                    
                    ! This ensures that the length is correct
                    distance_cyls(1:3) = cyls(i,6:8)-cyls(i,1:3)
                    length_cyls=norm2(distance_cyls)
                    unit_cyls=distance_cyls(1:3)/length_cyls
                    
                    cyls(i,6)= cyls(i,9)*unit_cyls(1)+cyls(i,1)
                    cyls(i,7)=cyls(i,9)*unit_cyls(2)+cyls(i,2)
                    cyls(i,8)=cyls(i,9)*unit_cyls(3)+cyls(i,3)
                    
                    write(50,*)i,cyls(i,1:8)
                    if(debug)write(*,*) 'cylinder coordinates: ',i,cyls(i,:)
                    if(debug)write(*,*) 'length check: ',cyls(i,9),norm2(cyls(i,6:8)-cyls(i,1:3))
                    
                endif
            endif
            !end do

            !add bubble around cylinder to check faster
            !center point:
            spheres(i,1:3)=((cyls(i,6:8)+cyls(i,1:3))/2.d0)
            !radius (half the length plus some)
            spheres(i,4)=(cyls(i,9)/2.d0)+0.5
            if(debug)write(*,*)'sphere bubble: ',spheres(i,:)

            do while((cyls(i,3) .gt. 0) .and. (cyls(i,8) .gt. 0) .and. (nn_count(i,1) .eq. 0) )
                if (nn_count(i,2) .ne. 0) then
                    write(*,*) 'Error in Drop/Contact'
                    stop
                end if
                if(debug)write(*,*)'checking contacts'
                if(debug)write(*,*) '[Main] Precontact length check: ',cyls(i,9),norm2(cyls(i,6:8)-cyls(i,1:3))
!                CALL contacts(i,j,k,cylinder_num,debug,spheres,cyls,nn_count,nn,r,d,p,q,contact,pr,qr,pq,buffer)
!!!Begin new contacts addition
                CALL new_contact(max_rad,spheres,cyls,i,cylinder_num,num_contact,contact_pt,curr_pt,i_pt,i_num)
                if(debug)write(*,*) '[Main] Postcontact length check: ',cyls(i,9),norm2(cyls(i,6:8)-cyls(i,1:3))
                nn_count(i,:)=0
                do cont_i=1,num_contact
                    nn_count(i,cont_i)=i_num(cont_i)
                end do
                r=contact_pt(1,:)
                
                    
                
!!!End
!                if (debug) then
!                    write(*,*)'Re-entered program cylinders place from contacts'
!                    write(*,*)'bug',i,cyls(i,3),nn_count(i,1),nn(i,1,:),nn(i,2,:)
!                end if
                !write(*,*) 'new_contact: i_num - ',i_num
                !write(*,*) '   contact 1: contact_pt,curr_pt,i_pt - ',contact_pt(1,:),curr_pt(1,:),i_pt(1,:)
                !write(*,*) '   contact 2: contact_pt,curr_pt,i_pt - ',contact_pt(2,:),curr_pt(2,:),i_pt(2,:)
                !write(*,*) '(old) contacts: nn_count(i,:) - ', nn_count(i,:)
                !write(*,*) '   r - ',r

                if (debug)write(*,*) '========================================================= i=',i
                do while ((cyls(i,3) .gt. 0) .and. (cyls(i,8) .gt. 0) .and. (nn_count(i,1) .ne. 0) .and. (nn_count(i,2) .eq. 0))
                    if(debug)write(*,*)'entering rotation'
                    !read (*,*)

                    !check slope
                    !find "rise"
                    if (cyls(i,3) .ge. cyls(i,8)) then
                        slopeone=cyls(i,3)-cyls(i,8)
                    else if (cyls(i,3) .lt. cyls(i,8)) then
                        slopeone=cyls(i,8)-cyls(i,3)
                    endif
                    !find run
                    slopetwo=sqrt((cyls(i,1)-cyls(i,6))**2+(cyls(i,2)-cyls(i,7))**2)
                    !Check for vertical
                    if (((slopeone/slopetwo) .gt. 15) .and. (nn_count(i,2) .eq. 0)) then
                        write(*,*)'rotated too far, dropping!',slopeone,slopetwo,slopeone/slopetwo
                        nn_count(i,1)=0
                    else !rotate!!
                        if(debug)write(*,*) '[Main] Pre-rotation length check: ',cyls(i,9),norm2(cyls(i,6:8)-cyls(i,1:3))
                        if (rotation_count(i,1).eq.5000)then
                            write(*,*)'CONTACT POINT 5000',r
                            write(*,*)'End Point1 5000',cyls(i,1:3)
                            write(*,*)'End Point2',cyls(i,6:8)
                            end if
                            if (rotation_count(i,1).gt.19997)then
                                write(*,*)'CONTACT POINT 19000 ',r
                                write(*,*)'End Point1',cyls(i,1:3)
                                write(*,*)'End Point2',cyls(i,6:8)
                            end if
                        CALL otter_cyl_rotate_wrapper(r,cyls,i,nn_count,i_pt(1,:),rotation_count,box,rotation_error)
                        if(debug)write(*,*) ' Out of Rotate: '
                        if (rotation_error)then
                            exit
                        end if
                        !read(*,*)

!                       CALL contacts(i,j,k,cylinder_num,debug,spheres,cyls,nn_count,nn,r,d,p,q,contact,pr,qr,pq,buffer)
!!!Begin new contacts addition
                        CALL new_contact(max_rad,spheres,cyls,i,cylinder_num,num_contact,contact_pt,curr_pt,i_pt,i_num)
                        !write(*,*) '   Out Contact '

                        nn_count(i,:)=0
                        do cont_i=1,num_contact
                            nn_count(i,cont_i)=i_num(cont_i)
                        end do
                        r=contact_pt(1,:)
!!!End
                        if (nn_count(i,1) .eq. 0) then
                            write(*,*) 'No longer in contact...!'
                            
                        end if

                        !Update spheres
                        spheres(i,1:3)=((cyls(i,6:8)+cyls(i,1:3))/2.d0)
                        !spheres(i,4)=(cyls(i,9)/2)+1

                    end if
                end do !rotate
                if (rotation_error)then
                    exit
                end if
                if ((nn_count(i,1) .eq. 0) .and. (cyls(i,3) .gt. 0) .and. (cyls(i,8) .gt. 0)) then
                    if (nn_count(i,2) .ne. 0) then
                        write(*,*) 'Error in Drop/Rotate/Contacts'
                        stop
                    end if
                    !drop
                    if(debug)write(*,*)'drop'
                    cyls(i,3)=cyls(i,3) - buffer
                    cyls(i,8)=cyls(i,8) - buffer
                    spheres(i,1:3)=((cyls(i,6:8)+cyls(i,1:3))/2.d0)
                    !spheres(i,4)=(cyls(i,9)/2)+1
                    if(debug)write(*,*) 'Cylinder Drop: ',cyls(i,3),cyls(i,8)
                endif
            enddo !drop loop?
            if (rotation_error)then
                write(*,*)'[Main] File error:',model ,'Press Enter to continue'
                read(*,*)
                exit
            end if
if (input) close(unit=16)
!Check if it's above the box.
            if ((cyls(i,3) .gt. box_range+max_rad) .and. (cyls(i,8) .gt. box_range+max_rad)) then
                i=i-1 ! remove current cyl
                cylinder_num=cylinder_num-9 ! decrement cylinder count
                if (debug)write(*,*) 'Removing cylinder above box...',cylinder_num
            else ! add cylinder and periodic images...
                periodic_added=0
                do periodic_i=1,3
                    do periodic_j=1,3
                        if (.not. ((periodic_i .eq. 2) .and. (periodic_j .eq. 2))) then
                            i=i+1
                            periodic_added=periodic_added+1
                            cyls(i,:)=cyls(i-periodic_added,:)
                            cyls(i,1)=cyls(i,1)+(periodic_i-2)*box_range
                            cyls(i,6)=cyls(i,6)+(periodic_i-2)*box_range
                            cyls(i,2)=cyls(i,2)+(periodic_j-2)*box_range
                            cyls(i,7)=cyls(i,7)+(periodic_j-2)*box_range
                        end if
                    end do
                end do
            end if

            i=i+1
            rotation_count(i,1)=0
            write(*,*) i
        enddo !sphere i until cylinder_num



        sum_rotation=0
        if (debug) write(*,*) 'Generating cylinders complete...'
        if(debug)then
            x=0
            do i=1,cylinder_num
                x=x+nn_count(i,1)
                sum_rotation=rotation_count(i,1)+sum_rotation
            enddo
            write(*,*)'nn count total',x
        endif
        average_rotation=sum_rotation/cylinder_num
        if(debug)write(*,*)'average rotation',average_rotation
203     format('_cone ',f0.3,',',f0.3,',',f0.3,' ',f0.3,' T ',f0.3,' A ',f0.3,',',f0.3,',',f0.3)
204     format('_box ',f0.3,',',f0.3,',',f0.3,' ',f0.3,',',f0.3,',',f0.3)
        write(out_unit,'(a)') '-osnap off'
        write(out_unit,'(a)') '_zoom -300,-300,-300 500,500,500'
        write(out_unit,'(a)') '_FACETRES 6'
        do i=1,cylinder_num!,9
            in_box=.true.
            !    do jj=1,3
            !      if ((cyls(i,jj) .lt. (0.0-cyls(i,4))) .or. (cyls(i,jj) .gt. (box+cyls(i,4)))) in_box=.false.
            !    end do
            !
            !    do jj=6,8
            !      if ((cyls(i,jj) .lt. (0.0-cyls(i,4))) .or. (cyls(i,jj) .gt. (box+cyls(i,4)))) in_box=.false.
            !    end do
            !
            ! cyls(i,4) = cyls(i,4) * 1.25
            ! cyls(i,5) = cyls(i,5) * 1.25
            !if(in_box) write(out_unit,203) cyls(i,1:8)
            if((in_box) .AND. (cyls(i,1) .GT. -15.0) .AND. (cyls(i,1) .LT. 45.0) .AND. (cyls(i,2) .LT. 45.0) &
            .AND. (cyls(i,2) .GT. -15.0) .AND. (cyls(i,3) .LT. 45.0) .AND. (cyls(i,3) .GT. -15.0) .AND. (cyls(i,6) .LT. 45.0) &
            .AND. (cyls(i,6) .GT. -15.0) .AND. (cyls(i,7) .LT. 45.0) .AND. (cyls(i,7) .GT. -15.0) &
            .AND. (cyls(i,8) .GT. -15.0) .AND. (cyls(i,8) .LT. 45.0)) write(out_unit,203) cyls(i,1:8)
        end do
        write(out_unit,'(a)') '_union all '

        !!  write(out_unit,'(a)') 'box 0,0,-5 100,100,5'
        !!  write(out_unit,'(a)') 'box 0,0,95 100,100,105'
        !  write(out_unit,'(a)') '_union ALL '
        !  write(out_unit,'(a)') 'FILEDIA 0'
        !  saveas='_saveas  '//trim(path)//trim(adjustl(full_name))//'.dwg'
        !  write(out_unit,'(a)') trim(adjustl(saveas))
        !  if ((lig_edge .eq. 'y') .or. (lig_edge .eq. 'Y')) then
            write(out_unit,'(a)') '_group create network  all '
            bigbox(1)=0-5000*box
            bigbox(2)=5000*box
            smallbox(1)=0-0.065*box
            smallbox(2)=1.065*box
        301 format('_box ',f0.4,',',f0.4,',',f0.4,' ',f0.4,',',f0.4,',',f0.4)
            write(out_unit,301) bigbox(1),bigbox(1),bigbox(1),bigbox(2),bigbox(2),bigbox(2)
            write(out_unit,'(a)') '_group create bigbox  last '
            if (small_box) then
                box_center=0.5*box
                box = 0.1*box
                box_bottom_corner=box_center-0.5*box
                box_top_corner=box_center+0.5*box
                write(out_unit,301) box_bottom_corner,box_bottom_corner,box_bottom_corner,&
                box_top_corner,box_top_corner,box_top_corner
            else
           
        !!    write(out_unit,301) 0.,0.,smallbox(1),box,box,smallbox(2)
            !write(out_unit,301) 0.,0.,0.,box,box,box
            write(out_unit,301) 2.,2.,2.,box,box,box
            end if !Small Box write
            write(out_unit,'(a)') '_group create smallbox  last '
            write(out_unit,'(a)') '_subtract g bigbox  g smallbox '
            write(out_unit,'(a)') '_subtract g network  g bigbox '
        !    saveas='_saveas  '//trim(path)//trim(adjustl(full_name))//'np.dwg'
        !    exportstl='_export '//trim(path)//trim(adjustl(full_name))//'np.stl all'
        !    write(out_unit,'(a)') '_move all  0,0,0 500,500,500'
        !    write(out_unit,'(a)') trim(saveas)
        !    write(out_unit,'(a,a)') trim(exportstl),' '

        !    end if
        !!if(debug .eqv. .false.) then
        !write(out_unit,'(a)') '_erase all '
        !write(out_unit,'(a)') '-purge g network y y'
        !write(out_unit,'(a)') '-purge g bigbox y y'
        !write(out_unit,'(a)') '-purge g smallbox y y'
        !!end if

        write(out_unit,'(a)') 'FILEDIA 1'
        close(unit=50)
        
        
    end do

end program





!Begin subroutines


SUBROUTINE nn_detailed(i,j,cylinder_num,k,pr,qr,r,cyls,nn_count,nn,debug)
!Purpose to save information of the detailed nn list
implicit none

!declare calling parameters
integer,intent(in) :: i,j,cylinder_num
real(kind=8),dimension(3,1) :: pr,qr,r
real(kind=8),dimension(cylinder_num,9),intent(in) :: cyls
real(kind=8),dimension(cylinder_num,1),intent(in) :: nn_count
logical :: debug
real(kind=8),dimension(cylinder_num,2,8), intent(out)::nn
real(kind=8),intent(in)::k


!declare local parameters
integer::onenn

!make sure values have increased
if (nn_count(i,1) .gt. k) then
!save data
if(debug) write(*,*) 'Detailed NN, count is > k'

if (nn_count(i,1) .eq. 1) then
  onenn=j
  nn(i,1,:) = j
  nn(i,1,1) = qr(1,1)-pr(1,1)
  nn(i,1,2) = qr(2,1)-pr(2,1)
  nn(i,1,3) = qr(3,1)-pr(3,1)
  nn(i,1,4) = r(1,1)
  nn(i,1,5) = r(2,1)
  nn(i,1,6) = r(3,1)
  nn(i,1,7) = cyls(j,3)
  nn(i,1,8) = cyls(j,8)
  if (debug) write(*,*) 'Detailed NN, count =1, particle: ', nn(i,1,1)
end if

if ((nn_count(i,1) .eq. 2) .and. (j .ne. onenn)) then
  nn(i,2,:) = j
  nn(i,2,1) = qr(1,1)-pr(1,1)
  nn(i,2,2) = qr(2,1)-pr(2,1)
  nn(i,2,3) = qr(3,1)-pr(3,1)
  nn(i,2,4) = r(1,1)
  nn(i,2,5) = r(2,1)
  nn(i,2,6) = r(3,1)
  nn(i,2,7) = cyls(j,3)
  nn(i,2,8) = cyls(j,8)

  if (debug) write(*,*) 'Detailed NN, count =2, particle: ', nn(i,1,1)
end if

if (debug) write (*,*) 'Detailed nnc list used: ', nn(i,:,:)
end if

RETURN
END SUBROUTINE nn_detailed

SUBROUTINE rotate(pq,i,j,cylinder_num,p,q,cyls,debug,contact,nn,nn_count,d,pr,qr,r,spheres,buffer)
    implicit none
    !declare calling parameters
    integer, intent(in)::i,j,cylinder_num
    real(kind=8),dimension(3,2),intent(inout)::p,q
    real(kind=8), dimension(cylinder_num,9), intent(inout)::cyls
    logical,intent(inout)::debug,contact
    real(kind=8),dimension(3,1),intent(inout)::pr,qr,r
    real(kind=8),dimension(cylinder_num,4),intent(inout)::spheres
    real(kind=8),dimension(cylinder_num,2),intent(inout)::nn_count
    real(kind=8),intent(inout)::d,buffer
    real(kind=8),dimension(cylinder_num,2,8), intent(out)::nn
    real(kind=8),dimension(3,3),intent(inout)::pq

    !declare dummy parameters
    real(kind=8)::pmag,qmag,theta,bmajor,a,phixy,phiyz,deltx,delty,deltz,delta,deltb,magx,magy,magz,dp,k
    real(kind=8),dimension(cylinder_num,3)::rcpb,rcpt
    real(kind=8)::arclength,pi=3.14159
    real(kind=8),dimension(3,1)::sp
    real(kind=8),dimension(3)::circle,cylinder
    logical::try
    real(kind=8),dimension(3)::slope,cyl_center
    real(kind=8)::x,y

    if (i .gt. 1) then
        bmajor=0.0
        try=.false.

        !calculate complementary angle between cylinder vectors
        pmag=sqrt(p(1,2)**2+p(2,2)**2+p(3,2)**2)
        !if (try)write(*,*)'pmag',pmag
        qmag=sqrt(q(1,2)**2+q(2,2)**2+q(3,2)**2)
        !if (try)write(*,*)'qmag',qmag
        theta=abs((pi/2)-acos((dot_product(p(1:3,2),q(1:3,2)))/pmag/qmag))
        if (try)write(*,*)'pmag',pmag,'qmag',qmag,'angle between cylinders in radians:',theta

        !major and minor axes:
        !if (try)write(*,*) 'J, radius',j,cyls(i,4)
        bmajor=cyls(j,4)
        !if (try)write(*,*)'major axes, b:',bmajor
        a=bmajor*cos(theta)
        if (try)write(*,*)'J, radius',j,cyls(j,4),'major axes, b:',bmajor,'minor axes, a:',a

!1. compute rcp,b,* for delta phi (rcpb: cp on bottom cylinder)
        rcpb(j,1:3)=r(1:3,1)
        if (try)write(*,*)'rcp,b:',rcpb(j,:)

    !calculate the change in the a and b axes after one rotation
        delta=a-(a*cos(.02))
        deltb=bmajor-(bmajor*cos(.02))
        if (try)write(*,*)'changes in a and b:',delta,deltb

    !find anlge of vector with respect to x (xy) and y (yz) axis
        phixy=atan(q(2,2)/q(1,2))
        phiyz=atan(q(3,2)/q(2,2))
        if (try)write(*,*)'phixy:',phixy
        if (try)write(*,*)'phiyz:',phiyz

    !find ratio of change in magnitude from displacement
        magx=(abs(q(1,2)))/(abs(q(1,2))+abs(q(2,2))+abs(q(3,2)))
        magy=(abs(q(2,2)))/(abs(q(1,2))+abs(q(2,2))+abs(q(3,2)))
        magz=(abs(q(3,2)))/(abs(q(1,2))+abs(q(2,2))+abs(q(3,2)))

        if (try)write(*,*)'magx:',magx,'magy:',magy,'magz:',magz

    !find changes in the cps
        deltx=magx*delta*cos(phixy)
        delty=magy*delta*sin(phixy)
        deltz=magz*delta*sin(phiyz)
        if (try)write(*,*)'delt-x,y,z:',deltx,delty,deltz

    !find new rcpb after rotation
        if (abs(r(1,1)-cyls(i,1)) .gt. abs(r(1,1)-cyls(i,6))) then
            rcpb(j,1)=rcpb(j,1)-deltx
        else if (abs(r(1,1)-cyls(i,1)) .le. abs(r(1,1)-cyls(i,6))) then
            rcpb(j,1)=rcpb(j,1)+deltx
        endif

        if (abs(r(2,1)-cyls(i,2)) .gt. abs(r(2,1)-cyls(i,7))) then
            rcpb(j,2)=rcpb(j,2)-delty
        else if (abs(r(2,1)-cyls(i,2)) .le. abs(r(2,1)-cyls(i,7))) then
            rcpb(j,2)=rcpb(j,2)+delty
        endif

        rcpb(j,3)=rcpb(j,3)-deltz-deltb
        if (try)write(*,*)'rcp,b,*:',rcpb(j,:)

!2.Compute arclength
        arclength=sqrt((a**2*(sin(.02))**2)+(bmajor**2*(cos(.02))**2))
        if (try)write(*,*)'arclength:',arclength

!3. Compute rcp,t,* for deltaphi (rcpt: cp on top)
        rcpt(i,1:3)=r(1:3,1)
        if (try)write(*,*)'rcp,t:',rcpt(i,:)

    !find anlge of vector with respect to x (xy) and y (yz) axis
        phixy=atan(p(2,2)/p(1,2))
        phiyz=atan(p(3,2)/p(2,2))
        if (try)write(*,*)'phixy:',phixy
        if (try)write(*,*)'phiyz:',phiyz

    !find ratio of change in magnitude from displacement
        magx=(abs(p(1,2)))/(abs(p(1,2))+abs(p(2,2))+abs(p(3,2)))
        magy=(abs(p(2,2)))/(abs(p(1,2))+abs(p(2,2))+abs(p(3,2)))
        magz=(abs(p(3,2)))/(abs(p(1,2))+abs(p(2,2))+abs(p(3,2)))

        if (try)write(*,*)'magx:',magx
        if (try)write(*,*)'magy:',magy
        if (try)write(*,*)'magz:',magz

        deltx=magx*arclength*cos(phixy)
        delty=magy*arclength*sin(phixy)
        deltz=magz*arclength*sin(phiyz)
        if (try)write(*,*)'delt-x,y,z:',deltx,delty,deltz

    !find new rcpt after rotation
        rcpt(i,1)=rcpt(i,1)+deltx
        rcpt(i,2)=rcpt(i,2)+delty
        rcpt(i,3)=rcpt(i,3)-deltz
        if (try)write(*,*)'rcp,t,*:',rcpt(i,:)

!4. Shift the cylinder the change in contact point
        if (abs(r(1,1)-cyls(i,1)) .gt. abs(r(1,1)-cyls(i,6))) then
            cyls(i,1)=cyls(i,1)-deltx
            cyls(i,6)=cyls(i,6)-deltx
            if(try)write(*,*)'shift x1',cyls(i,1),cyls(i,6)
        else if (abs(r(1,1)-cyls(i,1)) .le. abs(r(1,1)-cyls(i,6))) then
            cyls(i,1)=cyls(i,1)+deltx
            cyls(i,6)=cyls(i,6)+deltx
            if(try)write(*,*)'shift x2',cyls(i,1),cyls(i,6)
        endif

        if (abs(r(2,1)-cyls(i,2)) .gt. abs(r(2,1)-cyls(i,7))) then
            cyls(i,2)=cyls(i,2)-delty
            cyls(i,7)=cyls(i,7)-delty
            if(try)write(*,*)'shift y1',cyls(i,2),cyls(i,7)
        else if (abs(r(2,1)-cyls(i,2)) .le. abs(r(2,1)-cyls(i,7))) then
            cyls(i,2)=cyls(i,2)+delty
            cyls(i,7)=cyls(i,7)+delty
            if(try)write(*,*)'shift y2',cyls(i,2),cyls(i,7)
        endif

        cyls(i,3)=cyls(i,3)-deltz
        cyls(i,8)=cyls(i,8)-deltz
        if(try)write(*,*)'shift z',cyls(i,3),cyls(i,8)
        if (try)write(*,*)'cylinder position change:',cyls(i,:)


!5. make the cylinder tangent to the ellipse

!find the center of the bottom cylinder where the new contact point is
        cyl_center(1)=cyls(j,1)
        cyl_center(2)=cyls(j,2)
        cyl_center(3)=r(3,1)-cyls(j,4)

!find inverse of the slope of the center to the new cp
        slope(1:3)=(1/(rcpt(i,1:3)-cyl_center(1:3)))
        if(debug)write(*,*)'rotate, slope: ',slope

!find ratio of cylinder before(x) and after(y) the contact point
        x=(sqrt(abs(r(1,1)-cyls(i,1)))**2+(abs(r(2,1)-cyls(i,2)))**2+(abs(r(1,1)-cyls(i,1)))**2)/cyls(i,9)
        !if(debug)write(*,*)'rotate,x: ',x
        y=1-x
        if(debug)write(*,*)'rotate,x: ',x,'rotate,y: ',y

!find new starting coordinates
        cyls(i,1)=rcpt(i,1)+(slope(1)*x)
        cyls(i,2)=rcpt(i,2)+(slope(2)*x)
        cyls(i,3)=rcpt(i,3)+(slope(3)*x)

!find new ending coordinates
        cyls(i,6)=rcpt(i,1)+(slope(1)*y)
        cyls(i,7)=rcpt(i,2)+(slope(2)*y)
        cyls(i,8)=rcpt(i,3)+(slope(3)*y)

!CALL contacts(i,j,k,cylinder_num,debug,spheres,cyls,nn_count,nn,r,d,p,q,contact,pr,qr,pq,buffer)
!if (try) write(*,*) 'rotate SUBROUTINE, called in contacts'

    endif
    RETURN

end subroutine rotate
