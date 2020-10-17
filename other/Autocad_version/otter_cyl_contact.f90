module otter_cyl_contact

implicit none

contains

    subroutine new_contact(max_rad,spheres,cyls,curr_cyl,cylinder_num,num_contact,contact_pt,curr_pt,i_pt,i_num)
        implicit none

        !calling parameters
        integer,intent(in) :: cylinder_num ! total number of cylinders...
        integer,intent(in) :: curr_cyl ! cylinder number in cyls list of falling cylinder
        real(kind=8),dimension(cylinder_num,4),intent(in) :: spheres  ! Each cylinder has a sphere centered on the cylinder to allow a quick rough touching check
        real(kind=8),dimension(cylinder_num,9),intent(in) :: cyls ! This is the data on the cylinders themselves: cyls(i,1:3) -beginning; cyls(i,4:5)-beginning rad, end rad; cyls(i,6:8)-end point;cyls(i,9)-length
        integer,intent(out) :: num_contact ! number of contact points
        real(kind=8),dimension(2,3),intent(out) :: contact_pt, curr_pt, i_pt ! contact point, point on curr_cyl/i_cyl axis a radius away from contact: not the axis point of closest approach, because of sphere touching approximation
        integer,dimension(2),intent(out) :: i_num ! the up to two contacting cylinder numbers in cyls list
        real(kind=8),intent(in) :: max_rad ! max_radius

        !local parameters
        integer :: i,j,k,dir_j
        real(kind=8),dimension(4) :: distance
        logical :: debug=.false.,contact,check_contact
        real(kind=8) :: cyl_sphere_sep,curr_rad_v,i_rad_v,pos_ratio,curr_rad,i_rad,d,g,radii_ratio,rad_curr,rad_i
        integer :: curr_cyl_num,i_cyl_num,curr_alpha,i_alpha
        real(kind=8),dimension(3) :: curr_cyl_v,i_cyl_v,curr_center,i_center,curr_min,curr_max,i_min,i_max
        real(kind=8),dimension(3) :: pclose,qclose
        real(kind=8),dimension(3,2) :: p,q
        real(kind=8),dimension(3,3) :: pq
        integer::indx(2),n
        real(kind=8),dimension(2,2)::a
        real(kind=8),dimension(2,1)::b
        
        !initialize
        i=1
        num_contact=0
        do dir_j=1,3    ! find max/min extent of curr_cyl in 3 directions
            curr_min(dir_j)=cyls(curr_cyl,dir_j)
            curr_max(dir_j)=cyls(curr_cyl,dir_j+5)
            if (curr_min(dir_j) .gt. curr_max(dir_j)) then
                curr_min(dir_j)=cyls(curr_cyl,dir_j+5)
                curr_max(dir_j)=cyls(curr_cyl,dir_j)
            end if
        end do
        !check previous cylinders until 2 are in contact...
        do while ((i .lt. curr_cyl) .and. (num_contact .lt. 2))
            !Proximity check
            check_contact=.true.
            dir_j=1
            do while ((dir_j .lt. 4) .and. (check_contact)) !Cycle through 3 directions, if no overlap in any direction, then don't check contact....
                i_min(dir_j)=cyls(i,dir_j) ! same as above, finding max/min extent of i'th cylinder
                i_max(dir_j)=cyls(i,dir_j+5)
                if (i_min(dir_j) .gt. i_max(dir_j)) then
                    i_min(dir_j)=cyls(i,dir_j+5)
                    i_max(dir_j)=cyls(i,dir_j)
                end if
                !Checks if max_i is below min_curr, or min_i is above max_curr, in either case NO CONTACT, go on to next dir...
                if ((i_max(dir_j) .lt. curr_min(dir_j)-2.d0*max_rad) .or. (i_min(dir_j) .gt. &
                    & curr_max(dir_j)+2.d0*max_rad)) check_contact = .false.
                dir_j=dir_j+1
            end do

            if (check_contact) then !Do a detailed check only if all ranges overlap
                if(debug) write(*,*) 'Buffer spheres overlap, checking for contact: ',curr_cyl,i
                !initialize for line-to-line distance check:
                contact = .false.
                p(:,:) = 0.d0
                q(:,:)=0.d0
                pq(:,:)=0.d0

                !place starting point into p/q arrays:
                p(1:3,1)=cyls(curr_cyl,1:3) !p_1
                q(1:3,1)=cyls(i,1:3) !q_1
                !place direction vectors into p/q arrays:
                p(1:3,2)=cyls(curr_cyl,6:8)-cyls(curr_cyl,1:3) ! p(:,2)=p_1->p_2
                q(1:3,2)=cyls(i,6:8)-cyls(i,1:3) ! q(:,2) = q_1->q_2

                !solve for pq vector matrix:
                pq(1:3,1)=q(1:3,2)          ! pq(:,1) = q_1->q_2
                pq(1:3,2)=-p(1:3,2)         ! pq(:,2) = - (p_1->p_2)
                pq(1:3,3)=q(1:3,1)-p(1:3,1) ! pq(:,3) = p_1->q_1

                !if (debug) write(*,*) 'P:',p(:,:),'Q:',q(:,:),'PQ:',pq(:,:)

                !build equation matrix a:
                a(1,1)=dot_product(pq(1:3,1),p(1:3,2))  ! q_1->q_2 DOT -(p_1->p_2)
                a(1,2)=dot_product(pq(1:3,2),p(1:3,2))  ! -(p_1->p_2) DOT -(p_1->p_2)
                a(2,1)=dot_product(pq(1:3,1),q(1:3,2))  ! q_1->q_2 DOT q_1->q_2
                a(2,2)=dot_product(pq(1:3,2),q(1:3,2))  ! -(p_1->p_2) DOT q_1->q_2

                !build right hand solution matrix b:
                b(1,1)=-dot_product(pq(1:3,3),p(1:3,2)) ! p_1->q_1 DOT -(p_1->p_2)
                b(2,1)=-dot_product(pq(1:3,3),q(1:3,2)) ! p_1->q_1 DOT q_1->q_2

                !if (debug) write(*,*) 'Augmented matrix A generated', a
                !if (debug) write(*,*) 'Right hand solution matrix B generated', b

                !call lu decomposition subroutines
                CALL ludcmp(a,n,indx,g)
                CALL lubksb(a,n,indx,b)

                !plug in solution vector to PQ
                pq(1:3,1)=pq(1:3,1)*b(1,1)  ! q_1->q_2 SCALE BY s
                pq(1:3,2)=pq(1:3,2)*b(2,1)  ! -(p_1->p_2) SCALE BY t


                pclose(:)=p(:,1)-pq(:,2)    ! p1 - -(p_1->p_2)*t = p3
                qclose(:)=q(:,1)+pq(:,1)    ! q1 + q_1->q_2*s = q3

				if (debug) write(*,*)'[new_contact] pclose', pclose
				if (debug) write(*,*)'[new_contact] qclose', qclose

                !find pq vector by adding each row togehter (stored in first column)
                pq(1:3,1)=pq(1:3,1)+pq(1:3,2)+pq(1:3,3) ! q_1->q_2*s + -(p_1->p_2)*t + p_1->q_1 = p_3->q_3


                !find length of vector pq:
                d=sqrt((pq(1,1))**2+(pq(2,1))**2+(pq(3,1))**2)
                if (debug) write(*,*)i, '[new_contact] d',d

                !if (debug) write(*,*) 'distance between p and q:', d
                !if (debug) write(*,*) 'Cylinder radii:',i,cyls(i,4),j,cyls(j,4)

                !if length is less than the radii of cylinders then check to make sure points are on the line
                !Now with FRUSTRUMS
                if ((b(1,1) .le. 1.) .and. (b(1,1) .ge. 0.) .and. (b(2,1) .le. 1.) .and. (b(2,1) .ge. 0.)) then
                    rad_curr=cyls(curr_cyl,4)+( cyls(curr_cyl,5)-cyls(curr_cyl,4) )*b(2,1)
                    rad_i=cyls(i,4)+( cyls(i,5)-cyls(i,4) )*b(1,1)
                    if (debug) write(*,*) '[new_contact] rad_curr,rad_i: ',rad_curr,rad_i
                    if (d .le. (rad_curr+rad_i)) then
                        contact=.true.
                    else
                        !Check if ends touch - not implemented yet
                        if (debug) write(*,*) '[new_contact] contact is false'
                    end if

                    if (contact) then !cylinders overlap, update nearest neighbor information
                        num_contact=num_contact+1
                        if (debug) write(*,*) '[new_contact] contact is true',num_contact
                        radii_ratio=rad_curr/(rad_curr+rad_i) ! Relative radius of i (p) compared to j (q)
                        if (debug) write(*,*) '[new_contact] radii_ratio',radii_ratio
                        contact_pt(num_contact,:)=(1-radii_ratio)*pclose(:)+radii_ratio*qclose(:) ! Weighted contact point, not just midpoint
                        if (debug) write(*,*) '[new_contact] Contact point:',contact_pt(num_contact,:)

                        i_pt(num_contact,:)=qclose(:)
                        i_num(num_contact)=i
                        curr_pt(num_contact,:)=pclose(:)
                    end if

                end if
!~                  if (.not.contact.and.(d.le. (cyls(i,4)+cyls(i,4))))then
                 
!~                  call edge_contact(cyls,p,q,i,j,cylinder_num,qclose,contact)
!~                      if (debug) write(*,*) '[new_contact] Contact=', contact
                     
!~                  end if
!                !!ERROR HERE - DOES NOT ACCOUNT FOR FRUSTRUMS!
!                if (d .le. (cyls(curr_cyl,4)+cyls(i,4))) then
!                    !if (debug) write(*,*) 'Contact points P,Q:',pclose,qclose
!                    contact=.false.
!                    if ((b(1,1) .le. 1.) .and. (b(1,1) .ge. 0.) .and. (b(2,1) .le. 1.) .and. (b(2,1) .ge. 0.)) then
!                        contact=.true.
!                    else
!                        !!!!!NEED TO ADD IF ENDS TOUCH SIDES OR ENDS!!!
!                        if (debug) write(*,*) 'contact is false'
!                    endif
!
!                    if (contact) then !cylinders overlap, update nearest neighbor information
!                        !if (debug) write(*,*) 'contact is true'
!                        num_contact=num_contact+1
!                        !!ERROR HERE - DOES NOT ACCOUNT FOR FRUSTRUMS!
!                        radii_ratio=cyls(curr_cyl,4)/(cyls(curr_cyl,4)+cyls(i,4)) ! Relative radius of i (p) compared to j (q)
!                        contact_pt(num_contact,:)=(1-radii_ratio)*pclose(:)+radii_ratio*qclose(:) ! Weighted contact point, not just midpoint
!                        !if (debug) write(*,*) 'Contact point:',contact_pt(num_contact,:)
!
!                        i_pt(num_contact,:)=qclose(:)
!                        i_num(num_contact)=i
!                        curr_pt(num_contact,:)=pclose(:)
!                    end if
!                end if

            end if
            i=i+1
        end do

    end subroutine

    subroutine alt_contact(max_rad,spheres,cyls,curr_cyl,cylinder_num,num_contact,contact_pt,curr_pt,i_pt,i_num)
        implicit none

        !calling parameters
        integer,intent(in) :: cylinder_num ! total number of cylinders...
        integer,intent(in) :: curr_cyl ! cylinder number in cyls list of falling cylinder
        real(kind=8),dimension(cylinder_num,4),intent(in) :: spheres  ! Each cylinder has a sphere centered on the cylinder to allow a quick rough touching check
        real(kind=8),dimension(cylinder_num,9),intent(in) :: cyls ! This is the data on the cylinders themselves: cyls(i,1:3) -beginning; cyls(i,4:5)-beginning rad, end rad; cyls(i,6:8)-end point;cyls(i,9)-length
        integer,intent(out) :: num_contact ! number of contact points
        real(kind=8),dimension(2,3),intent(out) :: contact_pt, curr_pt, i_pt ! contact point, point on curr_cyl/i_cyl axis a radius away from contact: not the axis point of closest approach, because of sphere touching approximation
        integer,dimension(2),intent(out) :: i_num ! the up to two contacting cylinder numbers in cyls list
        real(kind=8),intent(in) :: max_rad ! max_radius

        !local parameters
        integer :: i,j,k,dir_j
        real(kind=8),dimension(4) :: distance
        logical :: debug=.false.,contact,check_contact
        real(kind=8) :: cyl_sphere_sep,curr_rad_v,i_rad_v,pos_ratio,curr_rad,i_rad
        integer :: curr_cyl_num,i_cyl_num,curr_alpha,i_alpha
        real(kind=8),dimension(3) :: curr_cyl_v,i_cyl_v,curr_center,i_center,curr_min,curr_max,i_min,i_max


        i=1
        num_contact=0
        do dir_j=1,3
            curr_min(dir_j)=cyls(curr_cyl,dir_j)
            curr_max(dir_j)=cyls(curr_cyl,dir_j+5)
            if (curr_min(dir_j) .gt. curr_max(dir_j)) then
                curr_min(dir_j)=cyls(curr_cyl,dir_j+5)
                curr_max(dir_j)=cyls(curr_cyl,dir_j)
            end if
        end do
        do while ((i .lt. curr_cyl) .and. (num_contact .lt. 2))
            !Proximity check
            check_contact=.true.
            dir_j=1
            do while ((dir_j .lt. 4) .and. (check_contact))
                i_min(dir_j)=cyls(i,dir_j)
                i_max(dir_j)=cyls(i,dir_j+5)
                if (i_min(dir_j) .gt. i_max(dir_j)) then
                    i_min(dir_j)=cyls(i,dir_j+5)
                    i_max(dir_j)=cyls(i,dir_j)
                end if
                if ((i_max(dir_j) .lt. curr_min(dir_j)-2.d0*max_rad) .or. (i_min(dir_j) .gt. &
                    & curr_max(dir_j)+2.d0*max_rad)) check_contact = .false.
                dir_j=dir_j+1
            end do

            if (check_contact) then
                if(debug) write(*,*) 'Buffer spheres overlap, checking for contact: ',curr_cyl,i

                !!!!! SET the distance between speres for cylinder contact check....
                cyl_sphere_sep=0.05

                curr_cyl_num=int((cyls(curr_cyl,9)-cyls(curr_cyl,4)-cyls(curr_cyl,5))/cyl_sphere_sep)+1
                i_cyl_num=int((cyls(i,9)-cyls(i,4)-cyls(i,5))/cyl_sphere_sep)+1
                curr_cyl_v=cyls(curr_cyl,6:8)-cyls(curr_cyl,1:3)
                i_cyl_v=cyls(i,6:8)-cyls(i,1:3)
                curr_rad_v=cyls(curr_cyl,5)-cyls(curr_cyl,4)
                i_rad_v=cyls(i,5)-cyls(i,4)
                curr_alpha=0
                contact=.false.
                !write(*,*) 'curr_cyl_v,i_cyl_v,curr_rad_v,i_rad_v: ',curr_cyl_v,i_cyl_v,curr_rad_v,i_rad_v
                do while ((curr_alpha .lt. curr_cyl_num) .and. (.not. contact))
                    pos_ratio=(real(curr_alpha,8))/(real(curr_cyl_num,8))
                    curr_center=cyls(curr_cyl,1:3)+(pos_ratio)*curr_cyl_v
                    curr_rad=cyls(curr_cyl,4)+(pos_ratio)*curr_rad_v
                    i_alpha=0
                    do while ((i_alpha .lt. i_cyl_num) .and. (.not. contact))
                        i_center=cyls(i,1:3)+(real(i_alpha,8)/real(i_cyl_num,8))*i_cyl_v
                        i_rad=cyls(i,4)+(real(i_alpha,8)/real(i_cyl_num,8))*i_rad_v
                        distance(1:3)=i_center-curr_center
                        distance(4)=norm2((/distance(1:3)/))
                        !write(*,*) ' dist,curr_rad,i_rad: ',distance,curr_rad,i_rad
                        !read(*,*)
                        if (distance(4)-curr_rad-i_rad .lt. 0.d0) contact=.true.
                        i_alpha=i_alpha+1
                    end do
                    curr_alpha=curr_alpha+1
                end do
                if (contact) then
                    if(debug) write(*,*) 'Contact!'
                    num_contact=num_contact+1
                    contact_pt(num_contact,:)=curr_center+distance(1:3)/distance(4)*(curr_rad-i_rad+distance(4))/2.d0
                    i_pt(num_contact,1:3)=i_center
                    i_num(num_contact)=i
                    curr_pt(num_contact,:)=curr_center
                end if
            end if
            i=i+1
        end do

    end subroutine

SUBROUTINE contacts(i,j,k,cylinder_num,debug,spheres,cyls,nn_count,nn,r,d,p,q,contact,pr,qr,pq,buffer)
    !Purpose: initially checks to see if any spheres overlap for first round
    implicit none
    !Declare callling parameters
    integer,intent(inout)::cylinder_num,i,j
    logical,intent(inout)::debug,contact
    real(kind=8),dimension(cylinder_num,4),intent(inout)::spheres
    real(kind=8),dimension(cylinder_num,9),intent(inout)::cyls
    integer,dimension(cylinder_num,2),intent(inout)::nn_count
    real(kind=8),dimension(cylinder_num,2,8),intent(inout)::nn
    real(kind=8),dimension(3) :: pclose,qclose
    real(kind=8),dimension(3,1),intent(out)::r,pr,qr
    real(kind=8),dimension(3,2),intent(inout)::p,q
    real(kind=8),intent(out)::d,k
    real(kind=8),intent(inout)::buffer
    real(kind=8),dimension(3,3),intent(out)::pq

    !Declare local parameters
    real(kind=8),dimension(4)::distance
    real(kind=8),dimension(2,2)::a
    real(kind=8),dimension(2,1)::b
    real(kind=8),dimension(3,2)::b_numeric
    integer::indx(2),n,l,h
    real(kind=8)::g,radii_ratio
    real(kind=8),dimension(3,1)::t,u
    logical::contactone,contacttwo

    debug=.false.

    !if(debug) write(*,*) 'contacts subroutine'
    k=nn_count(i,1)
    j=1
    if(debug) write(*,*) 'contacts subroutine, k: ',k

    do while ((j .lt. i) .and. (nn_count(i,2) .eq. 0))
        if (debug) write(*,*) 'Entered do loop in contacts subroutine',j,i

        !calc distances between sphere(i) and sphere(j)
        distance(1:3)=spheres(i,1:3)-spheres(j,1:3)
        if(debug) write(*,*)'spheres(in contact subroutine): ',spheres(i,:),spheres(j,:)
        if(debug) write(*,*)'distances: ',distance(1:3)
        distance(4)=(sqrt(distance(1)**2+distance(2)**2+distance(3)**2))-spheres(i,4)-spheres(j,4)
        if(debug) write(*,*) 'Contacts subroutine: NN loop: ',j,spheres(i,3),distance(4)

        !if the distance between them is negative then the spheres overlap and can check to see if cylinders overlap
        if (distance(4) .lt. 0) then
            if(debug) write(*,*) 'Checking against subsequent cylinder to see if it overlaps'

            !initialize arrays:
            contact = .false.
            p(:,:) = 0
            q(:,:)=0
            pq(:,:)=0

            !place starting point into p/q arrays:
            p(1:3,1)=cyls(i,1:3) !p_1
            q(1:3,1)=cyls(j,1:3) !q_1
            if (debug) write (*,*) 'initial points 1,2:',p(1:3,1),q(1:3,1)

            !place direction vectors into p/q arrays:
            p(1,2)=cyls(i,6)-cyls(i,1)
            p(2,2)=cyls(i,7)-cyls(i,2)
            p(3,2)=cyls(i,8)-cyls(i,3) ! p(:,2) = p_1->p_2

            q(1,2)=cyls(j,6)-cyls(j,1)
            q(2,2)=cyls(j,7)-cyls(j,2)
            q(3,2)=cyls(j,8)-cyls(j,3) ! q(:,2) = q_1->q_2

            !solve for pq vector matrix:
            pq(1:3,1)=q(1:3,2)          ! pq(:,1) = q_1->q_2
            pq(1:3,2)=-p(1:3,2)         ! pq(:,2) = - (p_1->p_2)
            pq(1:3,3)=q(1:3,1)-p(1:3,1) ! pq(:,3) = p_1->q_1

            if (debug) write(*,*) 'P:',p(:,:),'Q:',q(:,:),'PQ:',pq(:,:)

            !build equation matrix a:
            a(1,1)=dot_product(pq(1:3,1),p(1:3,2))  ! q_1->q_2 DOT -(p_1->p_2)
            a(1,2)=dot_product(pq(1:3,2),p(1:3,2))  ! -(p_1->p_2) DOT -(p_1->p_2)
            a(2,1)=dot_product(pq(1:3,1),q(1:3,2))  ! q_1->q_2 DOT q_1->q_2
            a(2,2)=dot_product(pq(1:3,2),q(1:3,2))  ! -(p_1->p_2) DOT q_1->q_2

            !build right hand solution matrix b:
            b(1,1)=-dot_product(pq(1:3,3),p(1:3,2)) ! p_1->q_1 DOT -(p_1->p_2)
            b(2,1)=-dot_product(pq(1:3,3),q(1:3,2)) ! p_1->q_1 DOT q_1->q_2

            if (debug) write(*,*) 'Augmented matrix A generated', a
            if (debug) write(*,*) 'Right hand solution matrix B generated', b

            !call lu decomposition subroutines
            CALL ludcmp(a,n,indx,g)
            CALL lubksb(a,n,indx,b)

            !plug in solution vector to PQ
            pq(1:3,1)=pq(1:3,1)*b(1,1)  ! q_1->q_2 SCALE BY s
            pq(1:3,2)=pq(1:3,2)*b(2,1)  ! -(p_1->p_2) SCALE BY t
            pclose(:)=p(:,1)-pq(:,2)    ! p1 - -(p_1->p_2)*t = p3
            qclose(:)=q(:,1)+pq(:,1)    ! q1 + q_1->q_2*s = q3

            !find pq vector by adding each row togehter (stored in first column)
            pq(1:3,1)=pq(1:3,1)+pq(1:3,2)+pq(1:3,3) ! q_1->q_2*s + -(p_1->p_2)*t + p_1->q_1 = p_3->q_3

            !find length of vector pq:
            d=sqrt((pq(1,1))**2+(pq(2,1))**2+(pq(3,1))**2)
            if (debug) write(*,*) 'distance between p and q:', d

            !if length is less than the radii of cylinders then check to make sure points are on the line
            if (debug) write(*,*) 'Cylinder radii:',i,cyls(i,4),j,cyls(j,4)
            if (d .le. (cyls(i,4)+cyls(j,4))) then
                if (debug) write(*,*) 'Contact point P:',pclose
                if (debug) write(*,*) 'Contact point Q:',qclose

                contact=.false.
                if ((b(1,1) .le. 1.) .and. (b(1,1) .ge. 0.) .and. (b(2,1) .le. 1.) .and. (b(2,1) .ge. 0.)) then
                    contact=.true.
                else
                    !!!!!NEED TO ADD IF ENDS TOUCH SIDES OR ENDS!!!
                    if (debug) write(*,*) 'contact is false'
                endif

                 !if they overlap find the contact point

!!!!! This is where the error is... this count increment occurred whether contact was true or not, that is
!!!!!   it occured every time a sphere overlapped, not just when they actually touched!  So I added the if statement!
                if (contact) then
                    if (debug) write(*,*) '[contact] contact is true'
                    radii_ratio=cyls(i,4)/(cyls(i,4)+cyls(j,4)) ! Relative radius of i (p) compared to j (q)
                    if(debug) write(*,*) '[contact] --------------------------------------------radii_ratio',radii_ratio
                    if(cylinder_num.eq.217) stop
                    r(1:3,1)=(1-radii_ratio)*pclose(:)+radii_ratio*qclose(:) ! Weighted contact point, not just midpoint
                    if (debug) write(*,*) '[contact] Contact point:',r
                    !cylinders overlap, update nearest neighbor information
                    if (nn_count(i,1) .eq. 0) then
                        nn_count(i,1) = j
                        if(debug) write (*,*) 'Contacts SUBROUTINE: NN_Count:   1', i, nn_count(i,1)
                    else if (nn_count(i,2) .eq. 0) then
                        nn_count(i,2) = j
                        if(debug) write (*,*) 'Contacts SUBROUTINE: NN_Count:   2', i, nn_count(i,:)
                    else
                        write(*,*) "Error in Contacts!"
                        stop
                    end if

                        !save information to detailed nn chart:

!!!! NOTE: As you have this written, k is sent to nn_detailed before being updated... is that right?
!!!! Should only happen when there is a contact, or every i,j pair?
                        !CALL nn_detailed(i,j,cylinder_num,k,pr,qr,r,cyls,nn_count,nn,debug)
                endif
                !    k=nn_count(i,1)
            endif ! check if d is less than sum of radii
             if (contact.and.(d.le. (cyls(i,4)+cyls(j,4))))then
                 call edge_contact(cyls,p,q,i,j,cylinder_num,qclose,contact)
             end if
        endif ! check if sphers overlap (dist(4) < 0)

        j=j+1
    enddo !j .lt. i

    if(debug)write(*,*)'End of contacts subroutine'
    RETURN

END SUBROUTINE contacts

SUBROUTINE ludcmp(a,n,indx,g)
!Purpose called to get reduced row echelon form
!implicit none
!declare calline paramter
integer::n,indx(n),NMAX,cylinder_num
real(kind=8)::g,a(2,2),TINY
!largest expected n, and a small number
parameter (NMAX=130, TINY=1.0e-20)

!declare local variables
integer::i,imax,j,k
real(kind=8)::aamax,dum,sum,vv(NMAX) !vv stores the implicit scaling of each row

n=2
g=1 !no row interchanges yet
do i=1,n  !loop over rows to get the implicit scaling information
    aamax=0.
    do j=1,n
      if(abs(a(i,j)) .gt. aamax) aamax=abs(a(i,j))
    enddo
    !if (aamax .eq. 0) pause 'singular matrix in ludcmp'  !No nonzero largest element
    vv(i) =1/aamax !save the scaling
enddo

do j=1,n  !this is the loop over columns of Crout's method
  do i=1,j-1 !covers up to, but not including i=j
        sum=a(i,j)
        do k=1,i-1
            sum=sum-a(i,k)*a(k,j)
        enddo
        a(i,j)=sum
    enddo
    aamax=0   !initialize for the search for largest pivot element
    do i=j,n  !covers i=j and i=j+1
         sum=a(i,j)
         do k=1,j-1
          sum=sum-a(i,k)*a(k,j)
         enddo
         a(i,j)=sum
         dum=vv(i)*abs(sum)  !figure of merit for the pivot
         if (dum .ge. aamax) then  !is it better than best so far?
          imax=i
          aamax=dum
      endif
    enddo

    if (j .ne. imax)then  !do we need to interchange rows?
      do k=1,n
          dum=a(imax,k)   !yes, do so...
          a(imax,k)=a(j,k)
          a(j,k)=dum
         enddo
         g=-g  !change parity of d
         vv(imax)=vv(j)   !also interchange scale factor
    endif
    indx(j)=imax
    if(a(j,j) .eq.0)a(j,j)=TINY

    if(j .ne. n) then  !divide by pivot element
        dum=1./a(j,j)
      do i=j+1,n
            a(i,j)=a(i,j)*dum
      enddo
    end if
enddo !go back for the next column in the reduction
RETURN
END SUBROUTINE ludcmp


SUBROUTINE lubksb(a,n,indx,b)
!purpose: uses backsubstitution to solve for s and t
!implicit none
!declare calling parameters
integer::n,indx(n)
real(kind=8)::a(2,2),b(n)

!declare local parameters
integer::i,ii,j,ll
real(kind=8)::sum

n=2
ii=0  !when ii is set to positive value, it will become the index of the first nonvanishing element of b
do i=1,n  ! do forward sub
    ll=indx(i)
    sum=b(ll)
    b(ll)=b(i)
    if (ii .ne. 0) then
      do j=ii,i-1
            sum=sum-a(i,j)*b(j)
      enddo
    else if (sum .ne. 0) then
        ii=i  !a nonzero element was encountered so have to do sums in the loop above
    endif
    b(i)=sum
enddo
do i=n,1,-1  !now do back sub
     sum=b(i)
     do j=i+1,n
      sum=sum-a(i,j)*b(j)
     enddo
     b(i)=sum/a(i,i)  !store a component to solution vector x
enddo
RETURN
END SUBROUTINE lubksb

SUBROUTINE edge_contact(cyls,p,q,i,j,cylinder_num,qclose,contact)
    !The purpose of this code is to determine if there is contact a small angles
    !It is possible that the normal contact routine misses these contacts
    !This code is implamented to catch these contacts

    !Global parameters
    integer,intent(in) :: cylinder_num ! total number of cylinders...
    real(kind=8),dimension(cylinder_num,9),intent(in) :: cyls
    real(kind=8),dimension(3,2),intent(in) ::p,q
    real(kind=8),dimension(3),intent(in) :: qclose 
    integer,intent(in) :: i,j

    !local parameters
    real(kind=8) :: phi,L2,dist
    real(kind=8),dimension(3) :: end_close
    logical  :: contact,debug
     
     debug=.true.
if(debug) write(*,*)'[edge_contact] i',i
if(debug) write(*,*)'[edge_contact] j',j
    !Angle between two fibers
    phi=acos(dot_product(p(:,2),q(:,2))/(norm2(p(:,2))*norm2(q(:,2))))
    if(debug) write(*,*)'[edge_contact] Closest approach: top fiber', p
    if(debug) write(*,*)'[edge_contact] Closest approach: bot fiber', q
    if(debug) write(*,*)'[edge_contact] Angle between two fibers', phi
    !Determine which end is closer to the
    if (norm2(cyls(i,1:3)-qclose).gt.norm2(cyls(i,6:8)-qclose))then
        end_close=cyls(i,6:8)
    else
        end_close=cyls(i,1:3)
    end if
    if(debug) write(*,*)'[edge_contact] Close end',end_close
    !Distance from close end to line of closest approach point 
    L2=norm2(end_close-qclose) 
    if(debug)write(*,*)'[edge_contact] L2', L2
!Distance from the close end point to the center line of the other fiber
    dist=L2*asin(phi)
    if (dist.lt.(2*cyls(i,4)))then
        contact=.true.
       if(debug) write(*,*)'Contact'
    end if
if(debug) write(*,*)'[edge_contact] dist', dist
if(debug) write(*,*)'[edge_contact] 2*cyls(i,4)',2*cyls(i,4)

end subroutine

end module

