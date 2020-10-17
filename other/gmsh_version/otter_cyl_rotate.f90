module otter_cyl_rotate

implicit none

contains

    subroutine otter_cyl_rotate_wrapper(r,cyls,i,nn_count,q_pt,rotation_count)
        implicit none

        !calling parameters
        real(kind=8),dimension(:,:),intent(inout)      :: cyls,rotation_count
        real(kind=8),dimension(3),intent(in)        :: r,q_pt
        integer,dimension(:,:),intent(in)   :: nn_count
        integer,intent(in)                             :: i
                !local parameters
        real(kind=8),dimension(3,3)     :: p,q
        real(kind=8),dimension(2,3)     :: rp
        real(kind=8),dimension(3)       :: p_uvec,qq
        real(kind=8),dimension(4)       :: center_of_gravity
        real(kind=8)                    :: pmag,qmag
        integer                 :: j,ii,z
        logical               :: debug=.false.
        
        rotation_count(i,1)=rotation_count(i,1)+1
        !if(rotation_count(i,1).lt.5000)z =1
        ! if (rotation_count(i,1).eq.5000.or.rotation_count(i,1).eq.10000)debug=.true.
        ! if (rotation_count(i,1).gt.19998)debug=.true.
        if (debug)write(*,*)'                                               '
        if (debug)write(*,*) '************** rotation_count***************',rotation_count(i,1)
        if (debug)write(*,*)'                                               '
        if (rotation_count(i,1).gt.50000)then
            write(*,*) 'ERROR: rotation count exceeds 30000'
            write(*,*)'[rotate]: z',z
            stop
        end if 
        center_of_gravity(1)=cyls(i,1)+(cyls(i,6)-cyls(i,1))/2
        center_of_gravity(2)=cyls(i,2)+(cyls(i,7)-cyls(i,2))/2
        center_of_gravity(3)=cyls(i,3)+(cyls(i,8)-cyls(i,3))/2
        center_of_gravity(4)=1.0 !only necessary for rotation math
        if (debug)write(*,*) '[rotate]:End Point One:', cyls(i,1:3)
        if (debug)write(*,*) '[rotate]:End Point Two:', cyls(i,6:8)
        if ((rotation_count(i,1).eq.5000).or.(rotation_count(i,1).eq.10000))then  
            z =z+1    
            write(*,*)'[rotate]: z',z
            cyls(i,1)=cyls(i,1)+(cyls(i,6)-cyls(i,1))
            cyls(i,2)=cyls(i,2)+(cyls(i,7)-cyls(i,1))
            cyls(i,3)=cyls(i,3)+(cyls(i,8)-cyls(i,3))+50
            cyls(i,6)=cyls(i,6)+(cyls(i,6)-cyls(i,1))
            cyls(i,7)=cyls(i,7)+(cyls(i,7)-cyls(i,2))
            cyls(i,8)=cyls(i,8)+(cyls(i,8)-cyls(i,3))+50
            !cyls(i,6:8)=cyls(i,6:8)+(cyls(i,6:8)-cyls(i,1:3))
            if (debug) write(*,*)'Shifting center of gravity'
            if (debug)write(*,*) 'End Point One:', cyls(i,1:3)
            if (debug)write(*,*) 'End Point Two:', cyls(i,6:8)
        endif 
        if ((abs(r(1)-center_of_gravity(1)).lt.0.9).and.(abs(r(2)-center_of_gravity(2)).lt.0.9)) then      
            cyls(i,1)=cyls(i,1)+0.1*(cyls(i,6)-cyls(i,1))
            cyls(i,2)=cyls(i,2)+0.1*(cyls(i,7)-cyls(i,1))
            cyls(i,3)=cyls(i,3)+0.1*(cyls(i,8)-cyls(i,3))
            cyls(i,6)=cyls(i,6)+0.1*(cyls(i,6)-cyls(i,1))
            cyls(i,7)=cyls(i,7)+0.1*(cyls(i,7)-cyls(i,2))
            cyls(i,8)=cyls(i,8)+0.1*(cyls(i,8)-cyls(i,3))
            if (debug) write(*,*)'Shifting center of gravity'
            if (debug)write(*,*) 'End Point One:', cyls(i,1:3)
            if (debug)write(*,*) 'End Point Two:', cyls(i,6:8)
        endif 
        !Set up v
        !values from main otter_cylinders
        !p and q : line "p" from p(1,1:3) to p(2,1:3), which is vector p(3,1:3)
        p(1,:)=cyls(i,1:3)
        p(2,:)=cyls(i,6:8)
        p(3,:)=p(2,:)-p(1,:)
        j=nn_count(i,1)
        q(1,:)=cyls(j,1:3)
        q(2,:)=cyls(j,6:8)
        q(3,:)=q(2,:)-q(1,:)
        pmag=norm2(p(3,:))
        qmag=norm2(q(3,:))
        p_uvec(:)=p(3,:)/pmag
        rp(1,:)=p(1,:)-r(:)
        rp(2,:)=p(2,:)-r(:)
        qq=dot_product((r-q(1,:)),q(3,:))/dot_product(q(3,:),q(3,:))*q(3,:)+q(1,:)
        if (debug)write(*,*) qq,q_pt
        center_of_gravity(1)=p(1,1)+(p(2,1)-p(1,1))/2
        center_of_gravity(2)=p(1,2)+(p(2,2)-p(1,2))/2
        center_of_gravity(3)=p(1,3)+(p(2,3)-p(1,3))/2
        center_of_gravity(4)=1.0 !only necessary for rotation math
        !read(*,*)
        
!        do ii=1,3
!            write(*,*) ' otter_cyl_rotate_wrapper: p(,',ii,') :  ',p(ii,:)
!            write(*,*) ' otter_cyl_rotate_wrapper: q(,',ii,') :  ',q(ii,:)
!            write(*,*) ' otter_cyl_rotate_wrapper: rp(,',ii,') : ',rp(ii,:)
!        end do
!        write(*,*) ' otter_cyl_rotate_wrapper: radii : ',cyls(i,4:5),cyls(j,4:5)

        CALL rotate_beck(r,qq,p,q,rp,p_uvec,pmag,qmag,center_of_gravity,debug)

!        if(debug) write(*,*) ' otter_cyl_rotate_wrapper: old cyl, new cyl: ', cyls(i,:),p(1,:),p(2,:)

        cyls(i,1:3)=p(1,:)
        cyls(i,6:8)=p(2,:)
        if (isnan(cyls(i,1))) stop 
    end subroutine

subroutine rotate_beck(r,qq,p,q,rp,p_uvec,pmag,qmag,center_of_gravity,debug)

    !calling parameters
    real(kind=8),dimension(3,3),intent(in)      :: q
    real(kind=8),dimension(3,3),intent(inout)   :: p
    real(kind=8),dimension(4),intent(inout)        :: center_of_gravity
    real(kind=8),dimension(2,3),intent(in)      :: rp
    real(kind=8),dimension(3),intent(in)        :: r,qq,p_uvec
    real(kind=8),intent(inout)                     :: pmag,qmag
    logical,intent(in)                     :: debug
    !local parameters
    real(kind=8),parameter              :: PI=4.0*atan(1.d0)
    
    logical                               :: align
    real(kind=8)                :: d,mag_unit_PQ,mag_unit_fiber,theta
    real(kind=8),dimension(4,4) :: T,T_inv,Rx, Rx_inv, Ry, Ry_inv,Rz,Rz_neg
    real(kind=8),dimension(3) :: unit_fiber,unit_PQ,CoG
    real(kind=8),dimension(4) :: cp,endpoint_one, endpoint_two,axis,qq_t
    
    qq_t(1:3)=qq(1:3)
    qq_t(4)=1.0
    CoG(1:3)=center_of_gravity(1:3)
    if (debug) write(*,*)'[Rotate_beck] CoG',CoG(1:3)
    if (debug) write(*,*)'[Rotate_beck] Center of gravity',center_of_gravity(1:3)
    ! Contact point to be translated
    cp(1)=r(1)
    cp(2)=r(2)
    cp(3)=r(3)
    cp(4)=1.0
    if (debug) write(*,*)'[Rotate_beck] contact point',cp(1:3)
!Defining X,Y,Z for end point 1
    endpoint_one(1)=p(1,1)
    endpoint_one(2)=p(1,2)
    endpoint_one(3)=p(1,3)
    endpoint_one(4)=1
    if (debug) write(*,*)'[Rotate_beck] endpoint_one',endpoint_one(1:3)
    !Defining X,Y,Z for end point 2
    endpoint_two(1)=p(2,1)
    endpoint_two(2)=p(2,2)
    endpoint_two(3)=p(2,3)
    endpoint_two(4)=1
    if (debug) write(*,*)'[Rotate_beck] endpoint_two',endpoint_two(1:3)
    !Define the translation and inverse translation 
    T=0.0
    T(1,1)=1.d0
    T(1,4)=-cp(1)
    T(2,4)=-cp(2)
    T(3,4)=-cp(3)
    T(2,2)=1.d0
    T(3,3)=1.d0
    T(4,4)=1.d0
    
    T_inv=0.0
    T_inv(1,1)=1.d0
    T_inv(1,4)=cp(1)
    T_inv(2,4)=cp(2)
    T_inv(3,4)=cp(3)
    T_inv(2,2)=1.d0
    T_inv(3,3)=1.d0
    T_inv(4,4)=1.d0
    
    !Translate end points
    endpoint_one=matmul(T,endpoint_one)
    endpoint_two=matmul(T,endpoint_two)
    center_of_gravity=matmul(T,center_of_gravity)
    qq_t=matmul(T,qq_t)
    cp=matmul(T,cp)
    if (debug) write(*,*)'Translated'
    if (debug) write(*,*)'[Rotate_beck] endpoint_one',endpoint_one(1:3)
    if (debug) write(*,*)'[Rotate_beck] endpoint_two',endpoint_two(1:3)
    if (debug) write(*,*)'[Rotate_beck] contact point',cp(1:3)
    !Defining a unit vector at the contact point pointing in the direction 
! of the top  fiber

unit_fiber(1)=endpoint_one(1)-endpoint_two(1)
unit_fiber(2)=endpoint_one(2)-endpoint_two(2)
unit_fiber(3)=endpoint_one(3)-endpoint_two(3)
mag_unit_fiber=sqrt(unit_fiber(1)**2+unit_fiber(2)**2+unit_fiber(3)**2)
unit_fiber=unit_fiber/mag_unit_fiber
if (debug) write(*,*)unit_fiber

!Defining a unit vector at the contact point pointing along PQ

unit_PQ(1:3)=qq_t(1:3)-cp(1:3)
mag_unit_PQ=sqrt(unit_PQ(1)**2+unit_PQ(2)**2+unit_PQ(3)**2)
unit_PQ=unit_PQ/mag_unit_PQ
if (debug) write(*,*)unit_PQ

!Calculating a unit vector at the contact point pointing along the axis of rotation

axis(1)= unit_fiber(2) *unit_PQ(3) - unit_fiber(3) *unit_PQ(2)
axis(2)= unit_fiber(3) *unit_PQ(1) - unit_fiber(1) *unit_PQ(3)
axis(3)= unit_fiber(1) *unit_PQ(2) - unit_fiber(2) *unit_PQ(1)
axis(4)=1.d0
if (debug) write(*,*)'Axis',axis
d=sqrt(axis(2)**2+axis(3)**2)

!Defining the rotation matrix
theta=0.01
Rx=0
Rx(1,1)=1
Rx(2,2)=axis(3)/d
Rx(2,3)=-axis(2)/d
Rx(3,2)=axis(2)/d
Rx(3,3)=axis(3)/d
Rx(4,4)=1

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

Rz=0
Rz(1,1)=cos(theta)
Rz(1,2)=-sin(theta)
Rz(2,1)=sin(theta)
Rz(2,2)=cos(theta)
Rz(3,3)=1
Rz(4,4)=1

Rz_neg=0
Rz_neg(1,1)=cos(-theta)
Rz_neg(1,2)=-sin(-theta)
Rz_neg(2,1)=sin(-theta)
Rz_neg(2,2)=cos(-theta)
Rz_neg(3,3)=1
Rz_neg(4,4)=1
align=.true.
!Rotation about X-axis
if ((axis(1).ne.1.0)) then
    endpoint_one=matmul(Rx,endpoint_one)
    endpoint_two=matmul(Rx,endpoint_two)
    center_of_gravity=matmul(Rx,center_of_gravity)
    cp=matmul(Rx,cp)
    axis=matmul(Rx,axis)
    align=.false.
    
    
    if (debug) write(*,*)'Rx'
    if (debug) write(*,*)'[Rotate_beck] endpoint_one',endpoint_one(1:3)
    if (debug) write(*,*)'[Rotate_beck] endpoint_two',endpoint_two(1:3)
    if (debug) write(*,*)'[Rotate_beck] contact point',cp(1:3)
        end if
        
endpoint_one=matmul(Ry,endpoint_one)
endpoint_two=matmul(Ry,endpoint_two)
center_of_gravity=matmul(Ry,center_of_gravity)
cp=matmul(Ry,cp)
axis=matmul(Ry,axis)
if (debug) write(*,*)'Ry'
if (debug) write(*,*)'[Rotate_beck] endpoint_one',endpoint_one(1:3)
if (debug) write(*,*)'[Rotate_beck] endpoint_two',endpoint_two(1:3)
if (debug) write(*,*)'[Rotate_beck] contact point',cp(1:3)
if (debug) write(*,*)'[Rotate_beck] Center of gravity',center_of_gravity(1:3)
if (debug) write(*,*)axis

   
endpoint_one=matmul(Rz,endpoint_one)
endpoint_two=matmul(Rz,endpoint_two)
center_of_gravity=matmul(Rz,center_of_gravity)
if (debug) write(*,*)'Rz'
if (debug) write(*,*)'[Rotate_beck] endpoint_one',endpoint_one(1:3)
if (debug) write(*,*)'[Rotate_beck] endpoint_two',endpoint_two(1:3)   


endpoint_one=matmul(Ry_inv,endpoint_one)
endpoint_two=matmul(Ry_inv,endpoint_two)
center_of_gravity=matmul(Ry_inv,center_of_gravity)
if (debug) write(*,*)'Ry inv'
if (debug) write(*,*)'[Rotate_beck] endpoint_one',endpoint_one(1:3)
if (debug) write(*,*)'[Rotate_beck] endpoint_two',endpoint_two(1:3)
if (debug) write(*,*)'[Rotate_beck] contact point',cp(1:3)

if (.not.align) then 
    endpoint_one=matmul(Rx_inv,endpoint_one)
    endpoint_two=matmul(Rx_inv,endpoint_two)
    center_of_gravity=matmul(Rx_inv,center_of_gravity)
        
    if (debug) write(*,*)'Rx inv'
    if (debug) write(*,*)'[Rotate_beck] endpoint_one',endpoint_one(1:3)
     if (debug) write(*,*)'[Rotate_beck] endpoint_two',endpoint_two(1:3)
end if

endpoint_one=matmul(T_inv,endpoint_one)
endpoint_two=matmul(T_inv,endpoint_two)
center_of_gravity=matmul(T_inv,center_of_gravity)
if (debug) write(*,*)'Translated Inv'
if (debug) write(*,*)'[Rotate_beck] endpoint_one',endpoint_one(1:3)
if (debug) write(*,*)'[Rotate_beck] endpoint_two',endpoint_two(1:3)

if (debug) write(*,*)'New'
if (debug) write(*,*)'[Rotate_beck] endpoint_one',endpoint_one(1:3)
if (debug) write(*,*)'[Rotate_beck] endpoint_two',endpoint_two(1:3)
if (debug) write(*,*)'[Rotate_beck] Center of gravity',center_of_gravity(1:3)
if(CoG(3).lt.center_of_gravity(3))then
    if (debug) write(*,*)'============================='
    if (debug) write(*,*)'=Rotated the wrong direction=' 
    if (debug) write(*,*)'============================='

    !Defining X,Y,Z for end point 1
    endpoint_one(1)=p(1,1)
    endpoint_one(2)=p(1,2)
    endpoint_one(3)=p(1,3)
    endpoint_one(4)=1
    if (debug) write(*,*)'[Rotate_beck] endpoint_one',endpoint_one(1:3)
    !Defining X,Y,Z for end point 2
    endpoint_two(1)=p(2,1)
    endpoint_two(2)=p(2,2)
    endpoint_two(3)=p(2,3)
    endpoint_two(4)=1
    if (debug) write(*,*)'[Rotate_beck] endpoint_two',endpoint_two(1:3)
    !Translate end points
    endpoint_one=matmul(T,endpoint_one)
    endpoint_two=matmul(T,endpoint_two)
    center_of_gravity=matmul(T,center_of_gravity)
    qq_t=matmul(T,qq_t)
    cp=matmul(T,cp)
    if (debug) write(*,*)'Translated'
    if (debug) write(*,*)'[Rotate_beck] endpoint_one',endpoint_one(1:3)
    if (debug) write(*,*)'[Rotate_beck] endpoint_two',endpoint_two(1:3)
    if (debug) write(*,*)'[Rotate_beck] contact point',cp(1:3)
    unit_fiber(1)=endpoint_one(1)-endpoint_two(1)
    unit_fiber(2)=endpoint_one(2)-endpoint_two(2)
    unit_fiber(3)=endpoint_one(3)-endpoint_two(3)
    mag_unit_fiber=sqrt(unit_fiber(1)**2+unit_fiber(2)**2+unit_fiber(3)**2)
    unit_fiber=unit_fiber/mag_unit_fiber
    if (debug) write(*,*)unit_fiber

    !Defining a unit vector at the contact point pointing along PQ

    unit_PQ(1:3)=qq_t(1:3)-cp(1:3)
    mag_unit_PQ=sqrt(unit_PQ(1)**2+unit_PQ(2)**2+unit_PQ(3)**2)
    unit_PQ=unit_PQ/mag_unit_PQ
    if (debug) write(*,*)unit_PQ

    !Calculating a unit vector at the contact point pointing along the axis of rotation

    axis(1)= unit_fiber(2) *unit_PQ(3) - unit_fiber(3) *unit_PQ(2)
    axis(2)= unit_fiber(3) *unit_PQ(1) - unit_fiber(1) *unit_PQ(3)
    axis(3)= unit_fiber(1) *unit_PQ(2) - unit_fiber(2) *unit_PQ(1)
    axis(4)=1.d0
    if (debug) write(*,*)'Axis',axis
    d=sqrt(axis(2)**2+axis(3)**2)
    align=.true.
    if ((axis(1).ne.1.0)) then
        endpoint_one=matmul(Rx,endpoint_one)
        endpoint_two=matmul(Rx,endpoint_two)
        center_of_gravity=matmul(Rx,center_of_gravity)
        cp=matmul(Rx,cp)
        axis=matmul(Rx,axis)
        
        align=.false.
        
        if (debug) write(*,*)'Rx'
        if (debug) write(*,*)'[Rotate_beck] endpoint_one',endpoint_one(1:3)
        if (debug) write(*,*)'[Rotate_beck] endpoint_two',endpoint_two(1:3)
        if (debug) write(*,*)'[Rotate_beck] contact point',cp(1:3)
         end if
    endpoint_one=matmul(Ry,endpoint_one)
    endpoint_two=matmul(Ry,endpoint_two)
    center_of_gravity=matmul(Ry,center_of_gravity)
    cp=matmul(Ry,cp)
    axis=matmul(Ry,axis)
    if (debug) write(*,*)'Ry'
    if (debug) write(*,*)'[Rotate_beck] endpoint_one',endpoint_one(1:3)
    if (debug) write(*,*)'[Rotate_beck] endpoint_two',endpoint_two(1:3)
    if (debug) write(*,*)'[Rotate_beck] contact point',cp(1:3)
    if (debug) write(*,*)'[Rotate_beck] Center of gravity',center_of_gravity(1:3)
    if (debug) write(*,*)axis
    if (debug) write(*,*)'[Rotate_beck] endpoint_two',endpoint_two(1:3)

    endpoint_one=matmul(Rz_neg,endpoint_one)
    endpoint_two=matmul(Rz_neg,endpoint_two)
    center_of_gravity=matmul(Rz_neg,center_of_gravity)
    if (debug) write(*,*)'Rz_neg'
     if (debug) write(*,*)'[Rotate_beck] endpoint_one',endpoint_one(1:3)
    if (debug) write(*,*)'[Rotate_beck] endpoint_two',endpoint_two(1:3)
        endpoint_one=matmul(Ry_inv,endpoint_one)
        endpoint_two=matmul(Ry_inv,endpoint_two)
        center_of_gravity=matmul(Ry_inv,center_of_gravity)
    if (debug) write(*,*)'Ry inv'
    if (debug) write(*,*)'[Rotate_beck] endpoint_one',endpoint_one(1:3)
    if (debug) write(*,*)'[Rotate_beck] endpoint_two',endpoint_two(1:3)
    if (debug) write(*,*)'[Rotate_beck] contact point',cp(1:3)
        
    if (.not.align) then
    endpoint_one=matmul(Rx_inv,endpoint_one)
    endpoint_two=matmul(Rx_inv,endpoint_two)
    center_of_gravity=matmul(Rx_inv,center_of_gravity)
        
    if (debug) write(*,*)'Rx inv'
    if (debug) write(*,*)'[Rotate_beck] endpoint_one',endpoint_one(1:3)
    if (debug) write(*,*)'[Rotate_beck] endpoint_two',endpoint_two(1:3)
    end if
        
    endpoint_one=matmul(T_inv,endpoint_one)
    endpoint_two=matmul(T_inv,endpoint_two)
    center_of_gravity=matmul(T_inv,center_of_gravity)
    if (debug) write(*,*)'Translated Inv'
    if (debug) write(*,*)'[Rotate_beck] endpoint_one',endpoint_one(1:3)
    if (debug) write(*,*)'[Rotate_beck] endpoint_two',endpoint_two(1:3)
        
    if (debug) write(*,*)'New'
    if (debug) write(*,*)'[Rotate_beck] endpoint_one',endpoint_one(1:3)
    if (debug) write(*,*)'[Rotate_beck] endpoint_two',endpoint_two(1:3)
endif


p(1,1:3)=endpoint_one(1:3)
p(2,1:3)=endpoint_two(1:3)


    end subroutine

end module
