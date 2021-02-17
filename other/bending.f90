module bending

     use otter_globals
     use otter_input
     use otter_fibers_globals
     use otter_math
     use bspline_module
     use bspline_kinds_module, only: wp, ip
 

    subroutine int_bending
    implicit none
    integer,parameter:: max_fibers=4
    integer(ip)::i,num_contacts,n,k
    real(wp),dimension(max_fibers,11,3)::fiber_path
    real(kind=8),dimension(max_fibers,8)::fibers
    real(kind=8),dimension(max_fibers,1)::bent
    real(kind=8),dimension(max_fibers,3)::p
    real(kind=8),dimension(max_fibers,11,1)::dist
    real(kind=8),allocatable::cp(:,:)
    bent=0.d0
   
    
    do i=1,4
    

    call fiber_mapping(fibers,i,fiber_path,max_fibers,num_contacts,cp,bent)
    if (i.gt.1)then
         do k=1,i-1
            do n=1,11
              call point_dist(fiber_path,fibers,i,k,n,max_fibers,dist)
            end do
            write(18,223)i,k,dist(k,1:11,1)
         end do
    end if
    deallocate(cp)
    end do
close(16)

    

subroutine fiber_mapping(fibers,i,fiber_path,max_fibers,num_contacts,cp,bent)
use bspline_module
use bspline_kinds_module, only: wp, ip
        implicit none
    !Inputs
    real(kind=8),dimension(4) :: p1,p2,endpoint1,endpoint2,cog
    integer,intent(in):: max_fibers,i,num_contacts
    real(kind=8),dimension(max_fibers,1)::bent
    real(kind=8),dimension(3,num_contacts)::cp
    !local parameters
    real(kind=8),parameter              :: PI=4.0*atan(1.d0)
    integer(ip)::num_pt,m
    logical                               ::debug,aligned
    real(kind=8)                :: d,mag_unit_fiber,length,length_original,lenght_tol
    real(kind=8),dimension(4,4) :: T,T_inv,Rx, Rx_inv, Ry, Ry_inv,Rz,Rz_neg
    real(kind=8),dimension(3) :: unit_fiber,unit_PQ,unit_cyl
    real(kind=8),dimension(4) :: axis,test,pt_p1_neg,pt_p1_pos,pt_p2_neg,pt_p2_pos,temp
    real(wp),dimension(5)::fx,y,dist
    real(wp),dimension(0:10) :: x_new,f_new_default
    real(wp),dimension(max_fibers,11,3),intent(out)::fiber_path
    real(kind=8),dimension(max_fibers,8)::fibers
    debug=.true.
    aligned=.false.
    num_pt=5
    
     endpoint1=1.d0
     endpoint2=1.d0
     endpoint1(1:3)=fibers(i,1:3)
     endpoint2(1:3)=fibers(i,6:8)
     p1=1.d0
     p1(1:3)=endpoint1(1:3)+((endpoint2(1:3)-endpoint1(1:3))/4.d0)
     p2=1.d0
     p2(1:3)=endpoint1(1:3)+((endpoint2(1:3)-endpoint1(1:3))/(4.0/3.0))
     cog=1.d0
     cog(1:3)=endpoint1(1:3)+((endpoint2(1:3)-endpoint1(1:3))/2.d0)
     dist=0.d0
     if(num_contacts.eq.2)then
         if (norm2(cp(1:3,1)-endpoint1(1:3)).lt.norm2(cp(1:3,2)-endpoint1(1:3))) then        
              dist(1)=norm2(cp(1:3,1)-endpoint1(1:3))
         else
              dist(1)=norm2(cp(1:3,2)-endpoint1(1:3))
         end if
         if (norm2(cog(1:3)-cp(1:3,1)).lt.norm2(cog(1:3)-cp(1:3,2))) then
              dist(3)=norm2(cog(1:3)-cp(1:3,1))
         else
              dist(3)=norm2(cog(1:3)-cp(1:3,2))
         end if
         if(norm2(endpoint2(1:3)-p1(1:3)).lt.norm2(endpoint2(1:3)-p2(1:3))) then
                 dist(5)=norm2(endpoint2(1:3)-p1(1:3))
         else
                 dist(5)=norm2(endpoint2(1:3)-p2(1:3))
         end if
         bent(i,1)=1
     end if
     if(num_contacts.eq.1)then 
             dist(1)=norm2(cp(1:3,1)-endpoint1(1:3))
             dist(2)=norm2(cp(1:3,1)-p1(1:3))
             dist(4)=norm2(cp(1:3,1)-p2(1:3))
             dist(5)=norm2(cp(1:3,1)-endpoint2(1:3))
             bent(i,1)=1
     end if
    if(debug)write(*,*)'length =',length
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


!Calculating a unit vector at the contact point pointing along the axis of rotation

axis(1)= unit_fiber(2)  
axis(2)= -unit_fiber(1) 
axis(3)= 0.d0
axis(4)=1.d0

axis(1:3)=(axis(1:3)/norm2(axis(1:3)))
!if (debug) if(debug)write(*,*)'Axis',axis
d=sqrt(axis(2)**2+axis(3)**2)
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

if (endpoint2(1).lt.0.0001)then
        write(*,*)'Aligned with y axis'
        aligned=.true.

end if
endpoint1(1:3)=endpoint1(1:3)+test(1:3)*dist(1)**2.d0
p1(1:3)=p1(1:3)+test(1:3)*dist(2)**2.d0
cog(1:3)=cog(1:3)+test(1:3)*dist(3)**2.d0
p2(1:3)=p2(1:3)+test(1:3)*dist(4)**2.d0
endpoint2(1:3)=endpoint2(1:3)+test(1:3)*dist(5)**2.d0
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
else

        fx(2)=endpoint1(1)
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
write(*,*)'>>>>>> X = ',fx
write(*,*)'y = ',y
        call sean_tests(fx,y,num_pt,x_new,f_new_default)
temp=1.d0
if (aligned)then
    do m=0,10
        write(*,*)'X=',x_new(m)
        write(*,*)'fx =',f_new_default(m)
        temp(2)=x_new(m)
        temp(1)=f_new_default(m)
        temp(3)=0.d0
        temp=matmul(Ry_inv,temp)
        if(d.gt.0.00001)temp=matmul(Rx_inv,temp)
        temp=matmul(T_inv,temp)
        fiber_path(i,m+1,1:3)=temp(1:3)
        write(*,*)'Fiber path = ',fiber_path(i,m+1,1:3)
    end do
else
        do m=0,10
           write(*,*)'X=',x_new(m)
           write(*,*)'fx =',f_new_default(m)
           temp(1)=x_new(m)
           temp(2)=f_new_default(m)
           temp(3)=0.d0
           temp=matmul(Ry_inv,temp)
           if(d.gt.0.00001)temp=matmul(Rx_inv,temp)
           temp=matmul(T_inv,temp)
           fiber_path(i,m+1,1:3)=temp(1:3)
           write(*,*)'Fiber path = ',fiber_path(i,m+1,1:3)
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
fibers(i,6:8)=endpoint2(1:3)

202 format('_test ',f0.3,',',f0.3,',',f0.3,' ',f0.3,',',f0.3,',',f0.3,' '&
           ,f0.3,',',f0.3,',',f0.3,' ',f0.3,',',f0.3,',',f0.3,' ',f0.3,','&
           ,f0.3,',',f0.3,' ',f0.3,',',f0.3,',',f0.3,' ',f0.3,',',f0.3,','&
           ,f0.3,' ',f0.3,',',f0.3,',',f0.3,' ',f0.3,',',f0.3,',',f0.3,' '&
           ,f0.3,',',f0.3,',',f0.3,' ',f0.3,',',f0.3,',',f0.3,' ',f0.3)
201 format('_five ',f0.3,',',f0.3,',',f0.3,' ',f0.3,',',f0.3,',',f0.3,' '&
            ,f0.3,',',f0.3,',',f0.3,' ',f0.3,',',f0.3,',',f0.3,' ',f0.3,','&
            ,f0.3,',',f0.3,' ',f0.3)
    write(16,201)endpoint1(1:3),p1(1:3),cog(1:3),p2(1:3),endpoint2(1:3),fibers(i,4)
    
end subroutine fiber_mapping

subroutine  sean_tests(fx,y,num_pt,x_new,f_new_default)

use bspline_module
use bspline_kinds_module, only: wp, ip

implicit none

    integer(ip),intent(in)::num_pt
    real(wp),dimension(num_pt),intent(in)::fx,y

    integer :: i    !! counter

     integer(ip),parameter :: kx  = 3    !! x bspline order
     integer(ip),parameter :: nx= 3      !! number of points in x dimension
     real(wp),dimension(num_pt):: x       !! [0,20,40,60,80,100]
     !real(wp),dimension(num_pt)::fx,y
     real(wp),dimension(num_pt)    :: fcn
     type(bspline_1d)          :: s_default
     real(wp),dimension(0:10),intent(out) :: x_new,f_new_default
     !,f_actual
     real(wp)                  :: xval,step_size
     integer(ip)               :: iflag
     x=fx
     fcn = y
     step_size=(x(num_pt)-x(1))/10

     !write(*,*)'x =',x
     !write(*,*)'y =',fcn
     call s_default%initialize(x,fcn,kx,iflag)  !default (not-a-knot)
     if (iflag/=0) error stop 'error initializing s_default'

     do i=0,10

        xval     = x(1)+i*step_size
        x_new(i) = xval

        !f_actual(i) = test_func(xval)

        call s_default%evaluate(xval,0_ip,f_new_default(i),iflag)
        if (iflag/=0) error stop 'error evaluating s_default'

     end do


end subroutine sean_tests
                                                               
subroutine point_dist(fiber_path,fibers,i,k,n,max_fibers,dist)
use bspline_kinds_module, only: wp, ip
implicit none 
    integer:: max_fibers
    real(wp),dimension(max_fibers,11,3)::fiber_path
    real(kind=8),dimension(3):: pt,endpoint1,endpoint2,AB,BE,AE,AB_cross_AE,AB_cross_BE
    real(kind=8),dimension(3)::p
    real(kind=8):: AB_dot_BE,AB_dot_AE,mag_AB,t
    real(kind=8),dimension(max_fibers,8)::fibers
    real(kind=8),dimension(max_fibers,11,1)::dist
    integer(ip)::i,n,k
    
    
    t=0
    pt(1)=fiber_path(k,n,1)
    pt(2)=fiber_path(k,n,2)
    pt(3)=fiber_path(k,n,3)
    endpoint1(1:3)=fibers(i,1:3)
    endpoint2(1:3)=fibers(i,6:8)
    write(*,*)pt

    AB=endpoint2-endpoint1
    BE=pt-endpoint2
    AE=pt-endpoint1
    write(*,*)AB
    write(*,*)BE
    write(*,*)AE
    AB_dot_BE=DOT_PRODUCT(AB,BE)
    AB_dot_AE=DOT_PRODUCT(AB,AE)
    if (AB_dot_AE.lt.0.d0) then
         dist(k,n,1)=sqrt(AE(1)**2.d0+AE(2)**2.d0+AE(3)**2d0)
    else
      if (AB_dot_BE.gt.0.d0)then
         dist(k,n,1)=sqrt(BE(1)**2.d0+BE(2)**2.d0+BE(3)**2d0)
      else
         AB_cross_AE(1)=AB(2)*AE(3)-(AB(3)*AE(2))
         AB_cross_AE(2)=AB(3)*AE(1)-(AB(1)*AE(3))
         AB_cross_AE(3)=AB(1)*AE(2)-(AB(2)*AE(1))
         write(*,*)'AB X AE = ',AB_cross_AE
         mag_AB=sqrt(AB(1)**2.d0+AB(2)**2.d0+AB(3)**2.d0)
         write(*,*)'mag.AB = ',mag_AB
         dist(k,n,1)=sqrt((AB_cross_AE(1)**2.d0+AB_cross_AE(2)**2.d0+AB_cross_AE(3)**2.d0))/mag_AB
         t=(pt(1)-endpoint1(1))/(AB(1))
         P=AB*t+endpoint1
      end if
    end if
    write(*,*)'AB . BE = ',AB_dot_BE
    write(*,*)'AB . AE = ',AB_dot_AE
    write(*,*)'Distance = ',dist(k,n,1)
    write(*,*)'T = ',t
    write(*,*)'P = ',p
end subroutine point_dist
                                                                    
end module bending