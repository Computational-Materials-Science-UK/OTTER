! Spheres - Skylars version...
!    made using beck's of version program OTTER
!    generates random structures based on spheres and particle overlap

program spheres_place

implicit none

! constants
real,allocatable,dimension(:,:) :: nn_count
real,allocatable,dimension(:,:) :: spheres, info
real,allocatable,dimension(:,:,:) :: nn
real,dimension(2) :: nn_avg
real(kind=8),dimension(4,4) :: t,rx,ry,rz,ti,rxi,ryi, newvec
integer :: num_models, sphere_num, model, i, j, jj,num_bad, out_unit=51, nn_c=47, nn_cd=46, magnitude, nn_num
character (len=18) :: out_name
character (len=200) :: path, cadpath, full_out_name, full_scr_name, full_nnc_name
character (len=200) ::full_nncd_name, saveas, exportsat, exportstl, full_name
real :: box, max_rad, min_rad, shift, sphere_range, box_range, sphere_buffer, drop,k, averagenn, overlapvolume, totalvolume
character :: lig_edge, cont
real,dimension(4) :: dist, distance
character (len=8),dimension(4) :: srf_char_array
real(kind=8),dimension(2) :: bigbox,smallbox
character (len=8) :: srf_char,num_char
logical :: good,debug,found,in_box
integer,dimension(8) :: values




! version history
integer,parameter :: version=3

! 1: random spheres with buffer, falland connect without rotating, nn count
! 2: attempt to rotate until touches 3 sheres
! 3: implentation of subroutines

debug=.false.

!batch information:

write (*,*) ' Welcome to Spheres v',version, '!'
write (*,*) ' '

out_unit=52 ! made for autoCAD

write (*,*) ' How many models do you want to make? '
read(*,*) num_models
write (*,*) ' What is the batch name of your models? '
read (*,*) out_name
write (*,*) ' What is the full path with trailing slash for your models? '
read(*,'(a)') path



!model information:

write(*,*) ' What is your desired box side length? '
read(*,*) box
write(*,*) ' Please enter minimum and maximum soft (outer) sphere radius: '
read(*,*) min_rad, max_rad
write(*,*) ' Please enter hard (inner) sphere volume: '
read(*,*) sphere_buffer
sphere_buffer = -sphere_buffer
if(debug) write(*,*) 'Sphere Buffer: ',sphere_buffer
write(*,*) ' Please enter the number of spheres: '
read(*,*) sphere_num
write(*,*) ' Ligamented edges? '
read (*,*) lig_edge

!ligamented edges:

shift=0.0


! Ligamented edges have ligaments that are truncted at the faces.  This is done by
! making a box with 2x the volume, and therefore cube-root-of-two-times the side length.
! To keep the density of spheres the same, the number of spheres is doubled.
! The "shift" places the corner of the unit box at the origin.
if ((lig_edge .eq. 'y') .or. (lig_edge .eq. 'Y')) then
  sphere_num=50*sphere_num
  box_range=3.0*box
  shift=0-1.0*box
else
  box_range=box
  shift=0.0
end if

!!structure:

write (*,*) ' Generating Structure(s)...'


! spheres:  sphere data for each of sphere_num spheres
!       x-pos, y-pos, z-pos, radius
allocate(spheres(sphere_num,4))

!  nn_count:  table showing how many spheres are overlapped with num_sphere
!   fields begin at 0 and increase by 1 each time a new sphere is added inside it's soft volume
!       sphere_num, number of nearest neighbors
allocate(nn_count(sphere_num,1))

!nn:  nearest neighbor table for each of sphere_num shperes
!       nn_num neighbors are stored with
!               atom_number, x-pos, y-pos, z-pos, distance, x-vec, y-vec, z-vec, overlap distance
!each atom gets 3 nn bc it falls until it touches 3, then all nn are only accounted for once
allocate(nn(sphere_num,3,8))

!info:  combines volume and nn data so it can be put into one sheet
!		stored with
!			model number, avgnn, total volume overlap
allocate(info(num_models,2))



!Full script
full_scr_name=trim(out_name)//".scr"
full_nnc_name=trim(out_name)//".txt"
open(unit=out_unit,file=trim(full_scr_name),status='new',action='write')
open(unit=47,file=trim(full_nnc_name),status='new',action='write')

!Begin building structures, one model at a time...e
do model=1, num_models

write(*,*) 'Model #',model,'...'
write(num_char,"(i4)") model !cast a four character version of the model number

full_name=trim(out_name)//trim(adjustl(num_char))
full_out_name=trim(out_name)//trim(adjustl(num_char))//".out"
full_nncd_name=trim(out_name)//trim(adjustl(num_char))//".plain"

open(unit=50,file=trim(full_out_name),status='new',action='write')
open(unit=46,file=trim(full_nncd_name),status='new',action='write')

!initialize arrays...
nn_count(:,1)=0.0
nn(:,:,:)=0.0

!get new random seed...
call date_and_time(values=values)
call srand(values(8))

sphere_range=max_rad-min_rad
if(debug) write(*,*) 'sphere_range: ',sphere_range

i=1
num_bad=0
do while (i .le. sphere_num)  !create spheres until sphere_num are created...

!randomly place a sphere above box
  !write(*,'(i3,a4)',ADVANCE='NO') i,'... '
  spheres(i,1)=rand()*box_range+shift
  spheres(i,2)=rand()*box_range+shift
  spheres(i,3)=rand()*box_range+shift+box
  spheres(i,4)=min_rad+rand()*sphere_range
  if(debug) write(*,*) 'Sphere loop: ',i,spheres(i,:)

do while ((spheres(i,3) .gt. 0) .and. (nn(i,3,1) .eq. 0.0))

CALL old_spheres (i,k,nn_count,spheres,nn,sphere_num,debug,sphere_buffer,overlapvolume)
if (debug) write(*,*) 'Re-entered program spheres_place from old_spheres'

    if (nn_count(i,1) .gt. k) then

      if ((nn(i,2,1) .eq. 0.0) .and. (nn(i,1,1) .gt. 0))then   !...if it only touches one sphere

            !rotate about sphere until it touches another sphere or both spheres are at the same height
            do while ((nn(i,2,1) .eq. 0.0) .and. (spheres(i,3) .gt. (nn(i,1,3))))

            CALL rotate_sphere (spheres,nn,sphere_num,i,debug,k,nn_count,sphere_buffer,overlapvolume)
            if (debug) write(*,*) 'Re-entered program spheres_place from rotate_sphere'


            end do  !rotate about sphere

            ! if the sphere needs to drop it no longer has a nearest neighbor
            if (spheres(i,3) .lt. distance(3)) then 
              nn(i,1,:) = 0
              nn_count(i,1) = nn_count(i,1) - 1
            end if

      else if ((nn(i,3,1) .eq. 0.0) .and. (nn(i,2,1) .gt. 0)) then !...if its touching two spheres

            !rotate about cone until it touches another sphere or they're at the same height 
            do while ((nn(i,3,1) .eq. 0.0) .and. (spheres(i,3) .gt. distance (3)))

            CALL rotate_cone (spheres,nn,sphere_num,i,debug,k,nn_count,sphere_buffer,overlapvolume)
            if (debug) write(*,*) 'Re-entered program spheres_place from rotate_cone'

            end do !rotate about cone

            ! if the sphere needs to drop it no longer has two nearest neighbors
            if (spheres(i,3) .lt. distance(3)) then
              nn(i,1,:) = 0
              nn(i,2,:) = 0
              nn_count(i,1) = nn_count(i,1) - 2
            end if

        end if !rotating

      else if ((nn_count(i,1) .eq. k) .and. (nn_count(i,1) .ne. 3)) then
        !drop 
        spheres(i,3)=spheres(i,3)+(rand()*sphere_buffer)
        if (debug) write (*,*) 'Sphere Drop: ',spheres(i,3)
    end if 

end do  ! sphere doesnt have 3nn and isnt touching the ground

i=i+1
end do  ! sphere i until sphere_num


if (debug) write(*,*) 'generating spheres comple'

write(50,*) ' Num spheres: ',sphere_num
write(50,*) ' Sphere buffer: ',sphere_buffer
write(50,*) ' Sphere rad min: ',min_rad
write(50,*) ' Sphere rad max: ',max_rad
write(50,*) ''
write(50,*) ' Sphere #, X, Y, Z, rad...'
do i=1,sphere_num
  write(50,*) i,spheres(i,:)
enddo

!write(47,*) ' Num spheres: ',sphere_num
!write(47,*) ''
!average nearest neighbor data
nn_avg = SUM (nn_count,DIM=1)
averagenn = nn_avg(1) / sphere_num
info(model,1) = averagenn
write(*,*) 'Average Number of Nearest Neighbors: ', averagenn

i=1
totalvolume = 0
!total volume overlap data
do while (i .le. sphere_num)
  j=1
  do while (j .le. 3)
    totalvolume = totalvolume + nn(i,j,8)
    if (debug) write(*,*) 'total volume do loop: ',i,j,nn(i,j,8)
    j=j+1
  end do
  i=i+1
end do
info(model,2) = totalvolume
write(*,*) ' Total volume of overlap in structure :', totalvolume


write(46,*) ' Num spheres: ',sphere_num
write(46,*) ''
write(46,*) 'Sphere #, Atom # of NN, x-pos, y-pos, z-pos, distance between spheres'
do i=1,sphere_num
     write(46,*) i,nn(i,1:3,1:4)
end do


srf_char_array(1)='    srf.'
srf_char_array(2)='   srf.0'
srf_char_array(3)='  srf.00'
srf_char_array(4)=' srf.000'

101 format('ic_surface sphere GEOM srf.',i5.5,' {{',3f14.10,'} {',3f14.10,'} ',f14.10,' 0 180}')
102 format('ic_surface cyl GEOM srf.',i5.5,' {{',3f14.10,'} {',3f14.10,'} ',2f14.10,' 1 1}')
201 format('_sphere ',a,' ',F14.1)

if(out_unit .eq. 52) then
!complex, parameter:: x = -175.
202 format('_sphere ',f0.3,',',f0.3,',',f0.3,' ',f0.3)
203 format('_cone ',f0.3,',',f0.3,',',f0.3,' ',f0.3,' T ',f0.3,' A ',f0.3,',',f0.3,',',f0.3)
204 format('_box ',f0.3,',',f0.3,',',f0.3,' ')
  write(out_unit,'(a)') '-osnap off'
  write(out_unit,'(a)') '_zoom -300,-300,-300 500,500,500'
  do i=1,sphere_num  
    in_box=.true.
    do jj=1,3
  if ((spheres(i,jj) .lt. (0.0-spheres(i,4))) .or. (spheres(i,jj) .gt. (box+spheres(i,4)))) in_box=.false.
    end do
  if (nn_count(i,1) .gt. 0) then
    if (in_box) write(out_unit,202) spheres(i,1:4)
  end if
  end do
  write(out_unit,'(a)') '_union all '
  

  write(out_unit,'(a)') 'FILEDIA 0'
  saveas='_saveas  '//trim(path)//trim(adjustl(full_name))//'.dwg'
  !write(out_unit,'(a)') trim(adjustl(saveas))
  if ((lig_edge .eq. 'y') .or. (lig_edge .eq. 'Y')) then
    write(out_unit,'(a)') '_group create network  all '
    bigbox(1)=-1*box
    bigbox(2)=3*box

301 format('_box ',f0.4,',',f0.4,',',f0.4,' ',f0.4,',',f0.4,',',f0.4)
    write(out_unit,301) bigbox(1),bigbox(1),bigbox(1),bigbox(2),bigbox(2),bigbox(2)
    write(out_unit,'(a)') '_group create bigbox  last '
    write(out_unit,301) 0.,0.,0.,box,box,box
    write(out_unit,'(a)') '_group create smallbox  last '
    write(out_unit,'(a)') '_subtract g bigbox  g smallbox '
    write(out_unit,'(a)') '_subtract g network  g bigbox '
    saveas='_saveas  '//trim(path)//trim(adjustl(full_name))//'np.dwg'
    exportstl='_export '//trim(path)//trim(adjustl(full_name))//'np.stl all'
    write(out_unit,'(a)') '_move all  0,0,0 500,500,500'
    write(out_unit,'(a)') trim(saveas)
    write(out_unit,'(a,a)') trim(exportstl),' '

  end if

  write(out_unit,'(a)') '_erase all '
  write(out_unit,'(a)') '-purge g network y y'
  write(out_unit,'(a)') '-purge g bigbox y y'
  write(out_unit,'(a)') '-purge g smallbox y y'



endif




close(unit=50)
close(unit=46)

end do  ! models

write(out_unit,'(a)') 'FILEDIA 1'


!write to .txt file with overall model information

write(47,*) ' Number of Structures :', num_models
write(47,*) ''
write(47,*) ' Model #,  Average # of NN, Total volume of overlap...'
do i=1, num_models
  write(47,*) i, info(i,:)
end do

close(unit=out_unit)
close(unit=47)

end program spheres_place

! subroutines

SUBROUTINE rotate_sphere (spheres,nn,sphere_num,i,debug,k,nn_count,sphere_buffer,overlapvolume)
!...Purpose: called in whenever a sphere is only touching one old sphere and needs to rotate
implicit none

! Declare calling parameters:
integer, intent(in) :: sphere_num,i
real, dimension(sphere_num,4), intent(inout) :: spheres
real, dimension(sphere_num,3,7), intent (inout) :: nn
real, dimension(sphere_num,1), intent(inout) :: nn_count
logical :: debug
real, intent (in) :: sphere_buffer
real, intent(inout) :: k
real, intent(out) :: overlapvolume

!Declare local variables
real :: x,y,z ! getting x,y location of second sphere
real :: a,b ! getting x,y location for first sphere
real :: co = 0.9996476952, si = 0.0174524064 ! rotation for 1 degree
real :: deltax, deltay ! info for changes in position
real :: d,e,f !transformations
real :: thetax,thetay ! for the roation matrix


if(debug) write(*,*) 'Sphere SUBROUTINE: ',i,spheres(i,:)

!starting position of the first sphere
a = nn(i,1,1)
b = nn(i,1,2)


if(debug) write(*,*) 'Sphere SUBROUTINE, first sphere position: ', a,b

!starting position of second sphere
x = spheres(i,1)
y = spheres(i,2)
z = spheres(i,3)


if(debug) write(*,*) ' Sphere SUBROUTINE, second sphere position: ', x,y


!information for the deltas
deltax = a-x
deltay = b-y

if(debug) write(*,*) 'Sphere SUBROUTINE, distance vectors: ',deltax,deltay

!transformations
if (abs(deltay) .ge. abs(deltax)) then 
  thetax = (deltay/(deltax+deltay))
  if (deltax .lt. 0) then
    CALL deltax_neg (x,y,z,thetax,d,e,f,debug)
    thetay = (deltax/(deltax+deltay))
    if (deltay .lt. 0) then
      CALL deltay_neg (thetay,d,e,f,debug)
    else if (deltay .gt. 0) then
      CALL deltay_pos (thetay,d,e,f,debug)
    end if

  else if (deltax .gt. 0) then
    CALL deltax_pos (x,y,z,thetax,d,e,f,debug)
    thetay = (deltax/(deltax+deltay))
    if (deltay .lt. 0) then
      CALL deltay_neg (thetay,d,e,f,debug)
    else if (deltay .gt. 0) then
      CALL deltay_pos (thetay,d,e,f,debug)
    end if
  end if

else if (abs(deltay) .lt. abs(deltax)) then 
  thetax = (deltax/(deltax+deltay))
  if (deltax .lt. 0) then
    CALL deltax_neg (x,y,z,thetax,d,e,f,debug)
    thetay = (deltax/(deltax+deltay))
    if (deltay .lt. 0) then
      CALL deltay_neg (thetay,d,e,f,debug)
    else if (deltay .gt. 0) then
      CALL deltay_pos (thetay,d,e,f,debug)
    end if

  else if (deltax .gt. 0) then
    CALL deltax_pos (x,y,z,thetax,d,e,f,debug)
    thetay = (deltax/(deltax+deltay))
    if (deltay .lt. 0) then
      CALL deltay_neg (thetay,d,e,f,debug)
    else if (deltay .gt. 0) then
      CALL deltay_pos (thetay,d,e,f,debug)
    end if
  end if

end if



if(debug) write(*,*) 'Sphere SUBROUTINE, matrix math: ', d,e,f

spheres(i,1) = d
spheres(i,2) = e
spheres(i,3) = f

if(debug) write (*,*) 'Sphere SUBROUTINE, final transformation: ',i,spheres(i,:)

CALL old_spheres (i,k,nn_count,spheres,nn,sphere_num,debug,sphere_buffer,overlapvolume)
if (debug) write(*,*) 'Sphere SUBROUTINE, called in old_spheres'

RETURN
END SUBROUTINE rotate_sphere


SUBROUTINE rotate_cone (spheres,nn,sphere_num,i,debug,k,nn_count,sphere_buffer,overlapvolume)
!...Purpose: called in whenever a sphere is only touching two old spheres and needs to rotate about a cone
implicit none

! Declare calling parameters:
integer, intent(in) :: sphere_num,i
real, dimension(sphere_num,4), intent(inout) :: spheres
real, dimension(sphere_num,3,7), intent (inout) :: nn
real, dimension(sphere_num,1), intent(inout) :: nn_count
logical :: debug
real, intent (in) :: sphere_buffer
real, intent(inout) :: k
real, intent(out) :: overlapvolume

!Declare local variables
real :: x,y,z ! getting location of previous sphere
real :: co = 0.9996476952, si = 0.0174524064 ! rotation for 1 degree
real :: a,b,c ! info for unit vector
real :: deltax, deltay ! info for changes in position
real :: d,e,f !transformations
real :: thetax,thetay ! for the roation matrix

if(debug) write(*,*) 'Cone SUBROUTINE: ',i,spheres(i,:)

!starting x,y,z of middle of cone vector
a = abs((nn(i,1,1)-nn(i,2,1))/2)
b = abs((nn(i,1,2)-nn(i,2,2))/2)
c = abs((nn(i,1,3)-nn(i,2,3))/2)

!starting position of second sphere
x = spheres(i,1)
y = spheres(i,2)
z = spheres(i,3)

!information for deltas
deltax = a-x
deltay = b-y

!transformations
if (abs(deltay) .ge. abs(deltax)) then 
  thetax = (deltay/(deltax+deltay))
  if (deltax .lt. 0) then
    CALL deltax_neg (x,y,z,thetax,d,e,f,debug)
    thetay = (deltax/(deltax+deltay))
    if (deltay .lt. 0) then
      CALL deltay_neg (thetay,d,e,f,debug)
    else if (deltay .gt. 0) then
      CALL deltay_pos (thetay,d,e,f,debug)
    end if

  else if (deltax .gt. 0) then
    CALL deltax_pos (x,y,z,thetax,d,e,f,debug)
    thetay = (deltax/(deltax+deltay))
    if (deltay .lt. 0) then
      CALL deltay_neg (thetay,d,e,f,debug)
    else if (deltay .gt. 0) then
      CALL deltay_pos (thetay,d,e,f,debug)
    end if
  end if

else if (abs(deltay) .lt. abs(deltax)) then 
  thetax = (deltax/(deltax+deltay))
  if (deltax .lt. 0) then
    CALL deltax_neg (x,y,z,thetax,d,e,f,debug)
    thetay = (deltax/(deltax+deltay))
    if (deltay .lt. 0) then
      CALL deltay_neg (thetay,d,e,f,debug)
    else if (deltay .gt. 0) then
      CALL deltay_pos (thetay,d,e,f,debug)
    end if

  else if (deltax .gt. 0) then
    CALL deltax_pos (x,y,z,thetax,d,e,f,debug)
    thetay = (deltax/(deltax+deltay))
    if (deltay .lt. 0) then
      CALL deltay_neg (thetay,d,e,f,debug)
    else if (deltay .gt. 0) then
      CALL deltay_pos (thetay,d,e,f,debug)
    end if
  end if

end if

!Take x,y,z out of parameters in deltay

if(debug) write(*,*) 'Cone SUBROUTINE, matrix math: ', d,e,f

spheres(i,1) = d
spheres(i,2) = e
spheres(i,3) = f

if(debug) write (*,*) 'Cone SUBROUTINE, final transformation: ',i,spheres(i,:)


CALL old_spheres (i,k,nn_count,spheres,nn,sphere_num,debug,sphere_buffer,overlapvolume)
if (debug) write(*,*) 'Cone SUBROUTINE, called in old_spheres'
RETURN
END SUBROUTINE rotate_cone

SUBROUTINE nn_detailed (nn_count,spheres,dist,i,j,k,nn,sphere_num,debug,overlapvolume)
!Purpose to save information of the detailed nn list
implicit none

!declare calling parameters
integer,intent(in) :: i,j,sphere_num
real, intent(in) :: k,overlapvolume
real,dimension(4) :: dist
real,dimension(sphere_num,4),intent(in) :: spheres
real,dimension(sphere_num,1),intent(in) :: nn_count
real,dimension(sphere_num,3,8), intent(out) :: nn
logical :: debug

!declare local parameters
!none
!make sure values have increased
if (nn_count(i,1) .gt. k) then
!save data
if (nn_count(i,1) .eq. 1) then
         nn(i,1,:) = j
         nn(i,1,1) = spheres(j,1)
         nn(i,1,2) = spheres(j,2)
         nn(i,1,3) = spheres(j,3)
         nn(i,1,4) = dist(4)
         nn(i,1,5) = dist(1)
         nn(i,1,6) = dist(2)
         nn(i,1,7) = dist(3)
         nn(i,1,8) = overlapvolume

         if (debug) write(*,*) 'Detailed NN: overlapvolume 1: ',nn(i,1,8)
         if (debug) write(*,*) 'Detailed NN: overlap distance: ', nn(i,1,4)
      else if (nn_count(i,1) .eq. 2) then
         nn(i,2,:) = j
         nn(i,2,1) = spheres(j,1)
         nn(i,2,2) = spheres(j,2)
         nn(i,2,3) = spheres(j,3)
         nn(i,2,4) = dist(4)
         nn(i,2,5) = dist(1)
         nn(i,2,6) = dist(2)
         nn(i,2,7) = dist(3)
         nn(i,2,8) = overlapvolume

         if (debug) write(*,*) 'Detailed NN: overlapvolume 1: ',nn(i,2,8)
         if (debug) write(*,*) 'Detailed NN: overlap distance: ', nn(i,2,4)
      else if (nn_count(i,1) .eq. 3) then
         nn(i,3,:) = j
         nn(i,3,1) = spheres(j,1)
         nn(i,3,2) = spheres(j,2)
         nn(i,3,3) = spheres(j,3)
         nn(i,3,4) = dist(4)
         nn(i,3,5) = dist(1)
         nn(i,3,6) = dist(2)
         nn(i,3,7) = dist(3)
         nn(i,3,8) = overlapvolume

         if (debug) write(*,*) 'Detailed NN: overlapvolume 1: ',nn(i,3,8)
         if (debug) write(*,*) 'Detailed NN: overlap distance: ', nn(i,3,4)
end if
if (debug) write (*,*) 'Detailed nnc list used: ', nn(i,:,:)

end if


RETURN
END SUBROUTINE nn_detailed

SUBROUTINE deltax_neg (x,y,z,thetax,d,e,f,debug)
!...Purpose: called in for transformation matrix, rotates about the y-axis ccw
implicit none

! Declare calling parameters:
real, intent (in) :: x,y,z
real, intent (in) :: thetax
real, intent (out) :: d,e,f
logical :: debug

if(debug) write (*,*) 'deltax_neg subroutine: x,y,z:', x,y,z
if(debug) write (*,*) 'deltax_neg subroutine: thetax:', thetax

d = x
e = (y*cos(thetax)-z*sin(thetax))
f = (y*sin(thetax)+z*cos(thetax))

if(debug) write (*,*) 'deltax_neg subroutine: d,e,f:', d,e,f

RETURN
END SUBROUTINE deltax_neg

SUBROUTINE deltax_pos (x,y,z,thetax,d,e,f,debug)
!...Purpose: called in for transformation matrix, rotates about the y-axis cw
implicit none

! Declare calling parameters:
real, intent (in) :: x,y,z
real, intent (in) :: thetax
real, intent (out) :: d,e,f
logical :: debug

if(debug) write (*,*) 'deltax_pos subroutine: x,y,z:', x,y,z
if(debug) write (*,*) 'deltax_pos subroutine: thetax:', thetax

d = x
e = (y*cos(thetax)+z*sin(thetax))
f = (-y*sin(thetax)+z*cos(thetax))

if(debug) write (*,*) 'deltax_pos subroutine: d,e,f:', d,e,f

RETURN
END SUBROUTINE deltax_pos

SUBROUTINE deltay_neg (thetay,d,e,f,debug)
!...Purpose: called in for transformation matrix, rotates about the x-axis ccw
implicit none

! Declare calling parameters:
real, intent (in) :: thetay
real, intent (inout) :: d,e,f
logical :: debug

if(debug) write (*,*) 'deltay_neg subroutine: x,y,z:', d,e,f
if(debug) write (*,*) 'deltay_neg subroutine: thetay:', thetay

d = (d*cos(thetay)-f*sin(thetay))
e = e
f = (e*sin(thetay)+f*cos(thetay))

if(debug) write (*,*) 'deltay_neg subroutine: d,e,f:', d,e,f

RETURN
END SUBROUTINE deltay_neg

SUBROUTINE deltay_pos (thetay,d,e,f,debug)
!...Purpose: called in for transformation matrix, rotates about the x-axis ccw
implicit none

! Declare calling parameters:
real, intent (in) :: thetay
real, intent (inout) :: d,e,f
logical :: debug

if(debug) write (*,*) 'deltay_pos subroutine: x,y,z:', d,e,f
if(debug) write (*,*) 'deltay_pos subroutine: thetay:', thetay

d = (d*cos(thetay)+f*sin(thetay))
e = e
f = (-e*sin(thetay)+f*cos(thetay))

if(debug) write (*,*) 'deltay_pos subroutine: d,e,f:', d,e,f

RETURN
END SUBROUTINE deltay_pos

SUBROUTINE old_spheres (i,k,nn_count,spheres,nn,sphere_num,debug,sphere_buffer,overlapvolume)
!..Purpose: called in to check new sphere against old sphere
implicit none

!Declare calling parameters
integer, intent(in) :: sphere_num,i
real, dimension(sphere_num,4), intent(in) :: spheres
real, dimension(sphere_num,3,8), intent (inout) :: nn
real, dimension(sphere_num,1), intent(inout) :: nn_count
logical :: debug
real, intent (in) :: sphere_buffer
real, intent (out) :: k, overlapvolume

!Declare local variables
real, dimension(4) :: dist
integer :: j

if (debug) write(*,*) 'Old Sphere SUBROUTINE'

k=nn_count(i,1)
j=1
do while (j .lt. i)
   !calc distances between sphere (i) and sphere (j)
   dist(1:3)=spheres(i,1:3)-spheres(j,1:3)
   dist(4)=sqrt(dist(1)**2+dist(2)**2+dist(3)**2)-spheres(i,4)-spheres(j,4)
        if(debug) write(*,*) 'Old_spheres SUBROUTINE: NN loop: ',j,dist(4)

   
    !Check and see if sphere being considered (i) is inside soft radius of test sphere (j)

    if((dist(4) .lt. 0) .and. (dist(4) .gt. sphere_buffer)) then
    !if it's inside the sphere and outside hard volume then nn count increases by one
      nn_count(i,1)=nn_count(i,1)+1
      nn_count(j,1)=nn_count(j,1)+1
      if(debug) write (*,*) 'Old_spheres SUBROUTINE: NN_Count:   ', i, nn_count(i,:)

      !save information to detailed nn chart
       CALL nn_detailed (nn_count,spheres,dist,i,j,k,nn,sphere_num,debug,overlapvolume)

      !Get volume overlapped
       CALL overlap_volume (i,j,k,nn,spheres,sphere_num,debug,overlapvolume,nn_count)

       k=nn_count(i,1)
    end if

j=j+1
end do ! j .lt. i

if (debug) write(*,*) 'End of Old Sphere SUBROUTINE'
RETURN
END SUBROUTINE old_spheres

SUBROUTINE overlap_volume (i,j,k,nn,spheres,sphere_num,debug,overlapvolume,nn_count)
!...Purpose: called to caluclate overlap volume between two spheres
implicit none

!Declare calling parameters
integer, intent(in) :: sphere_num,i,j
real, dimension(sphere_num,4), intent(in) :: spheres
real, dimension(sphere_num,3,8), intent (inout) :: nn
real, dimension(sphere_num,1), intent(in) :: nn_count
logical :: debug
real, intent (out) :: overlapvolume
real, intent(in) :: k

!Declare local variables
real :: volumeone, volumetwo, pi=3.14159, h

if (debug) write(*,*) 'Overlap: value of k: ',k
if (nn_count(i,1) .gt. k) then

if (nn_count(i,1) .eq. 1) then
  if (debug) write(*,*) 'Overlap: distance overlap 1: ',nn(i,1,4)
  h = (-0.5)*nn(i,1,4)
  if (debug) write(*,*) 'Overlap: h: ',h

else if (nn_count(i,1) .eq. 2) then
  if (debug) write(*,*) 'Overlap: distance overlap 2: ',nn(i,2,4)
  h = (-0.5)*nn(i,2,4)
  if (debug) write(*,*) 'Overlap: h: ',h

else if (nn_count(i,1) .eq. 3) then
  if (debug) write(*,*) 'Overlap: distance overlap 3: ',nn(i,3,4)
  h = (-0.5)*nn(i,3,4)
  if (debug) write(*,*) 'Overlap: h: ',h

end if

!Dome math of both spheres
volumeone = (.3333333)*pi*(h**2)*((3*spheres(i,4))-h)
if (debug) write(*,*) 'Overlap: volume of sphere 1: ',volumeone

volumetwo = (.3333333)*pi*(h**2)*((3*spheres(j,4))-h)
if (debug) write(*,*) 'Overlap: volume of sphere 2: ',volumetwo

overlapvolume = volumeone+volumetwo

!save overlap to detailed nn
if (nn_count(i,1) .eq. 1) then
  nn(i,1,8) = overlapvolume
  if (debug) write(*,*) 'Overlap: total overlap 1: ',nn(i,1,8)

else if (nn_count(i,1) .eq. 2) then
  nn(i,2,8) = overlapvolume
  if (debug) write(*,*) 'Overlap: total overlap 2: ',nn(i,2,8)

else if (nn_count(i,1) .eq. 3) then
  nn(i,3,8) = overlapvolume
  if (debug) write(*,*) 'Overlap: total overlap3: ',nn(i,3,8)

end if

end if
RETURN
END SUBROUtINE overlap_volume
