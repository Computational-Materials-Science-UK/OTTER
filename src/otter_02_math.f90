!!! OTTER - otter_math.f90
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
! otter_math - module providing basic math handling to the OTTER codebase
!	vb.1 - Working version, see acknowledgements in line.

! External Dependencies:
!	None
! Internal Dependencies:
!   otter_globals.mod

module otter_math

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Dependencies
    use otter_globals

    implicit none

    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Tools provided:
    !
    !> FUNCTION cross3(a, b) result(cross) - cross product of vectors in 3-space
    !> FUNCTION pt2line (x1, x2, x) result(d) - distance from a point (x) to a line through x1 and x2, all x in n-space, must be vectors of same length
    !> FUNCTION pt2seg(x1, x2, x) result(d) - distance from a point (x) to a *segment* (x1->x2), all x in n-space, must be vectors of same length
    !> recursive FUNCTION ax2om_d(a) result(res) - converts axis/angle to rotation matrix; a is four doubles, the first three the direction cosines (miller indices), and the fourth the angle to rotate in radians.

    FUNCTION cross3(a, b) result(cross)
        implicit none
        real(kind=DBL), DIMENSION(3) :: cross
        real(kind=DBL), DIMENSION(3), INTENT(IN) :: a, b
      
        cross(1) = a(2) * b(3) - a(3) * b(2)
        cross(2) = a(3) * b(1) - a(1) * b(3)
        cross(3) = a(1) * b(2) - a(2) * b(1)
    END FUNCTION cross3

    FUNCTION dotn(a,b) result(dot)
        implicit none
        real(kind=DBL) :: dot
        real(kind=DBL), DIMENSION(3), INTENT(IN) :: a, b

        dot=sum(a(:)*b(:))
    END FUNCTION dotn

    SUBROUTINE nearpt_pt2line3 (p,a,b,npt,dist)
        ! See https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Vector_formulation
        implicit none 
        real(kind=DBL),dimension(3),intent(IN) :: p,a,b
        real(kind=DBL),dimension(3),intent(OUT) :: npt 
        real(kind=DBL),intent(OUT) :: dist

        real(kind=DBL),dimension(3) :: n

        n=(b-a)/norm2((b-a))
        npt=(a-p)-dotn((a-p),n)*n
        dist=norm2(npt)

    END SUBROUTINE nearpt_pt2line3

    !!! Provided by BECK

    function cross2(a,b) result(cross)
        implicit None

        real(kind=DBL),dimension(2),intent(in) :: a,b
        real(kind=DBL) :: cross

        cross=a(1)*b(2)-a(2)*b(1)

    end function

    !*****************************************************************************************
    !> distance_from_point_to_line
    !> author: Jacob Williams
    !  date:8/2012
        !
        !  Compute the distance between the point X and the line defined
        !  by the two points X1 and X2.
        !
        !# References
        !  1. http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
    !
    function pt2line (x1, x2, x) result(d)

        implicit none

        real(kind=DBL),dimension(3),intent(in) :: x1
        real(kind=DBL),dimension(3),intent(in) :: x2
        real(kind=DBL),dimension(3),intent(in) :: x
        real(kind=DBL)                         :: d

        real(kind=DBL),dimension(3) :: x21
        real(kind=DBL) :: x21_mag

        x21 = x2 - x1
        x21_mag = norm2(x21)

        if (x21_mag/=0.d0) then
            d = norm2( cross3( x21, x1 - x ) ) / x21_mag
        else
            d = norm2(x1 - x)
        end if

    end function pt2line

    !*****************************************************************************************
    !> distance_from_point_to_line_segment
    !> author: Jacob Williams
    !  date:8/2012
        !
        !  Compute the distance between a line segment and a point.
        !
        !# References
        !  1. http://forums.codeguru.com/showthread.php?194400-Distance-between-point-and-line-segment
        !
    !@note x,x1,x2 should all be the same length

    function pt2seg(x1, x2, x) result(d)

        implicit none

        real(kind=DBL),dimension(:),intent(in) :: x1
        real(kind=DBL),dimension(:),intent(in) :: x2
        real(kind=DBL),dimension(:),intent(in) :: x
        real(kind=DBL)                         :: d

        real(kind=DBL),dimension(size(x1)) :: x12
        real(kind=DBL) :: s
        real(kind=DBL) :: x12_mag

        x12 = x1 - x2
        x12_mag = norm2(x12)

        if (x12_mag==0.d0) then
            d = norm2(x1 - x)
        else
            s = dot_product( x1-x, x12 ) / dot_product( x12, x12 )

            !if the projection is not on the segment,
            ! then use the closest end point:
            if (s<0.d0) then
                s = 0.d0
            else if (s>1.d0) then
                s = 1.d0
            end if

            d = norm2( x - x1 + s*x12 )
        end if

    end function pt2seg

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! The following subroutines were taken from rotations.f90
    ! https://github.com/marcdegraef/3Drotations/blob/master/src/f90/rotations.f90
    
    ! ###################################################################
    ! Copyright (c) 2013-2014, Marc De Graef/Carnegie Mellon University
    ! All rights reserved.
        !
        ! Redistribution and use in source and binary forms, with or without modification, are 
        ! permitted provided that the following conditions are met:
        !
        !     - Redistributions of source code must retain the above copyright notice, this list 
        !        of conditions and the following disclaimer.
        !     - Redistributions in binary form must reproduce the above copyright notice, this 
        !        list of conditions and the following disclaimer in the documentation and/or 
        !        other materials provided with the distribution.
        !     - Neither the names of Marc De Graef, Carnegie Mellon University nor the names 
        !        of its contributors may be used to endorse or promote products derived from 
        !        this software without specific prior written permission.
        !
        ! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
        ! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
        ! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
        ! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
        ! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
        ! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
        ! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
        ! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
        ! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
        ! USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
        ! ###################################################################

        !--------------------------------------------------------------------------
        ! EMsoft:rotations.f90
        !--------------------------------------------------------------------------
        !
        ! MODULE: rotations
        !
        !> @author Marc De Graef, Carnegie Mellon University
        !
        !> @brief everything that has to do with rotations and conversions between rotations
        !
        !> @details This file relies a lot on the relations listed in the book "Orientations
        !> and Rotations" by Adam Morawiec [Springer 2004].  I've tried to implement every
        !> available representation for rotations in a way that makes it easy to convert 
        !> between any pair.  Needless to say, this needs extensive testing and debugging...
        !>
        !> Instead of converting all the time between representations, I've opted to 
        !> "waste" a little more memory and time and provide the option to precompute all the representations.
        !> This way all representations are available via a single data structure.
        !>
        !> Obviously, the individual conversion routines also exist and can be called either in
        !> single or in double precision (using a function interface for each call, so that only
        !> one function name is used).  The conversion routines use the following format for their
        !> call name:  ab2cd, where (ab and cd are two-characters strings selected from the following
        !> possibilities: [the number in parenthesis lists the number of entries that need to be provided] 
        !>
        !> eu : euler angle representation (3)
        !> om : orientation matrix representation (3x3)
        !> ax : axis angle representation (4)
        !> ro : Rodrigues vector representation (3)
        !> qu : unit quaternion representation (4)
        !> ho : homochoric representation (3)
        !> cu : cubochoric representation (3).
        !>
        !> hence, conversion from homochoric to euler angle is called as ho2eu(); the argument of 
        !> each routine must have the correct number of dimensions and entries.
        !> All 42 conversion routines exist in both single and double precision.
        !>
        !> Some routines were modified in July 2014, to simplify the paths in case the direct conversion
        !> routine does not exist.  Given the complexity of the cubochoric transformations, all routines
        !> going to and from this representation will require at least one and sometimes two or three
        !> intermediate representations.  cu2eu and qu2cu currently represent the longest computation 
        !> paths with three intermediate steps each.
        !>
        !> In August 2014, all routines were modified to account for active vs. passive rotations,
        !> after some inconsistencies were discovered that could be traced back to that distinction.
        !> The default is for a rotation to be passive, and only those transformation rules have been
        !> implemented.  For active rotations, the user needs to explicitly take action in the calling 
        !> program.
        !>
        !> Testing: the program rotationtest.f90 was generated by an IDL script and contains all possible
        !> pairwise and triplet transformations, using a series of input angle combinations; for now, these
        !> are potentially problematic Euler combinations.
        !>
        !> The conventions of this module are:
        !>
        !> - all reference frames are right-handed and orthonormal (except for the Bravais frames)
        !> - a rotation angle is positive for a counterclockwise rotation when viewing along the positive rotation axis towards the origin
        !> - all rotations are interpreted in the passive way
        !> - Euler angles follow the Bunge convention, with phi1 in [0,2pi], Phi in [0,pi], and phi2 in [0,2pi]
        !> - rotation angles (in axis-angle derived representations) are limited to the range [0,pi]
        !> 
        !> To make things easier for the user, this module provides a routine to create a rotation 
        !> representation starting from an axis, described by a unit axis vector, and a rotation angle.
        !> This routine properly takes the sign of epsijk into account, and always produces a passive rotation.
        !> The user must explicitly take action to interpret a rotation a being active.
        !>
        !> @date 08/04/13 MDG 1.0 original
        !> @date 07/08/14 MDG 2.0 modifications to several routines (mostly simplifications)
        !> @date 08/08/14 MDG 3.0 added active/passive handling (all routines passive)
        !> @date 08/11/14 MDG 3.1 modified Rodrigues vector to 4 components (n and length) to accomodate Infinity
        !> @date 08/18/14 MDG 3.2 added RotateVector, RotateTensor2 routines with active/passive switch
        !> @date 08/20/14 MDG 3.3 completed extensive testing of epsijk<0 mode; all tests passed for the first time !
        !> @date 08/21/14 MDG 3.4 minor correction in om2ax to get things to work for epsijk>0 mode; all tests passed!
        !> @date 09/30/14 MDG 3.5 added routines to make rotation definitions easier
        !> @date 09/30/14 MDG 3.6 added strict range checking routines for all representations (tested on 10/1/14)
        !> @date 03/11/15 MDG 3.7 removed the RotVec_q routine, since it is surpassed by the new quat_Lp and quat_Lpstar routines
        !> @date 03/12/15 MDG 3.8 correction of Rodrigues representation for identity rotation -> [0,0,epsijk,0]
        !> @date 03/16/15 MDG 3.9 added quat_average routine
        !> @date 04/17/15 MDG 3.91 simplification to qu2eu routines
    !-------------------------------------------------------------------------- 

    !--------------------------------------------------------------------------
    ! Function: ax2om_d
    !
    !> @author Marc De Graef, Carnegie Mellon University
    !> @brief Axis angle pair to orientation matrix (double precision)
    !> @note verified 8/5/13.
    !> @param a axis angle pair (double precision)  
    !>  
    !> @date 8/04/13   MDG 1.0 original
    !--------------------------------------------------------------------------
    recursive function ax2om_d(a) result(res)

        IMPLICIT NONE

        real(kind=dbl),INTENT(IN)       :: a(4)         !< axis angle pair
        real(kind=dbl)                  :: res(3,3)
        
        real(kind=dbl)                  :: q, c, s, omc
        integer                         :: i

        c = dcos(a(4))
        s = dsin(a(4))
        omc = 1.D0-c

        do i=1,3
            res(i,i) = a(i)**2*omc + c
        end do

        q = omc*a(1)*a(2)
        res(1,2) = q + s*a(3)
        res(2,1) = q - s*a(3)

        q = omc*a(2)*a(3)
        res(2,3) = q + s*a(1)
        res(3,2) = q - s*a(1)

        q = omc*a(3)*a(1)
        res(3,1) = q + s*a(2)
        res(1,3) = q - s*a(2)

        ! See discussion in documentation about rotations.f90
        !if (epsijkd.eq.1.D0) res = transpose(res)

    end function ax2om_d

    !*****************************************************************************80
    !
    !! SEGMENTS_DIST_3D computes the distance between two line segments in 3D.
    !
        !  Discussion:
        !
        !    NOTE: The special cases for identical and parallel lines have not been
        !    worked out yet; those cases are exceptional, and so this code
        !    is made available in a slightly unfinished form!
        !
        !    A line segment is the finite portion of a line that lies between
        !    two points P1 and P2.
        !
        !    Given two line segments, consider the underlying lines on which
        !    they lie.
        !   
        !    A) If the lines are identical, then the distance between the line segments
        !    is 0, if the segments overlap, or otherwise is attained by the
        !    minimum of the distances between each endpoint and the opposing
        !    line segment.
        !
        !    B) If the lines are parallel, then the distance is either the distance
        !    between the lines, if the projection of one line segment onto
        !    the other overlaps, or otherwise is attained by the
        !    minimum of the distances between each endpoint and the opposing
        !    line segment.
        !
        !    C) If the lines are not identical, and not parallel, then there are 
        !    unique points PN and QN which are the closest pair of points on the lines.
        !    If PN is interior to [P1,P2] and QN is interior to [Q1,Q2],
        !    then the distance between the two line segments is the distance
        !    between PN and QN.  Otherwise, the nearest distance can be computed
        !    by taking the minimum of the distance from each endpoing to the
        !    opposing line segment.
        !
        !    Therefore, our computation first checks whether the lines are
        !    identical, parallel, or other, and checks for the special case
        !    where the minimum occurs in the interior.
        !
        !    If that case is ruled out, it computes and returns the minimum of:
        !
        !      Distance ( P1, [Q1,Q2] );
        !      Distance ( P2, [Q1,Q2] );
        !      Distance ( Q1, [P1,P2] );
        !      Distance ( Q2, [P1,P2] );
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license. 
        !
        !  Modified:
        !
        !    12 August 2006
        !
        !  Author:
        !
        !    John Burkardt
        !

    subroutine segments_dist_3d ( p1, p2, q1, q2, dist, pn, qn, sn, tn )

        !  Parameters:
        !
        !    Input, real(kind=DBL) P1(3), P2(3), the endpoints of the first
        !       segment.
        !    Input, real(kind=DBL) Q1(3), Q2(3), the endpoints of the second
        !       segment.
        !    Output, real(kind=DBL) DIST, the distance between the line segments.
        implicit none
        
        integer, parameter                 :: dim_num = 3
        integer                            :: close, i
        real(kind=DBL)                     :: a,b,c,d,e
        real(kind=DBL)                     :: det
        real(kind=DBL), intent(out)        :: dist
        real(kind=DBL), dimension(dim_num),intent(in)  :: p1,p2,q1,q2
        real(kind=DBL), dimension(dim_num),intent(out) ::   pn,qn
        real(kind=DBL), dimension(dim_num) :: u,v,w0
        real(kind=DBL), intent(out)        :: sn, tn
        real(kind=DBL), dimension(4)       :: un, dist2
        real(kind=DBL), dimension(4,3)     :: rn

        !  The lines are identical.
        !  THIS CASE NOT SET UP YET
        !
            ! if ( lines_exp_equal_3d ( p1, p2, q1, q2 ) ) then
            ! end if
            !
            !  The lines are not identical, but parallel
            !  THIS CASE NOT SET UP YET.
            !
            ! if ( lines_exp_parallel_3d ( p1, p2, q1, q2 ) ) then
            ! end if
            !
            !  C: The lines are not identical, not parallel.
        
        !  Let U = (P2-P1) and V = (Q2-Q1) be the direction vectors on
        !  the two lines.
        u(1:dim_num) = p2(1:dim_num) - p1(1:dim_num)
        v(1:dim_num) = q2(1:dim_num) - q1(1:dim_num)
    
        !> Let SN be the unknown coordinate of the nearest point PN on line 1,
        !  so that PN = P(SN) = P1 + SN * (P2-P1).
        !> Let TN be the unknown coordinate of the nearest point QN on line 2,
        !  so that QN = Q(TN) = Q1 + TN * (Q2-Q1).
        !> Let W0 = (P1-Q1).
        w0(1:dim_num) = p1(1:dim_num) - q1(1:dim_num)

        !  The vector direction WC = P(SN) - Q(TC) is unique (among directions)
        !  perpendicular to both U and V, so
        !    U dot WC = 0
        !    V dot WC = 0
        !  or, equivalently:
        !    U dot ( P1 + SN * (P2 - P1) - Q1 - TN * (Q2 - Q1) ) = 0
        !    V dot ( P1 + SN * (P2 - P1) - Q1 - TN * (Q2 - Q1) ) = 0
        !  or, equivalently:
        !    (u dot u ) * sn - (u dot v ) tc = -u * w0
        !    (v dot u ) * sn - (v dot v ) tc = -v * w0
        !  or, equivalently:
        !   ( a  -b ) * ( sn ) = ( -d )
        !   ( b  -c )   ( tc )   ( -e )
        a = dot_product ( u, u )
        b = dot_product ( u, v )
        c = dot_product ( v, v )
        d = dot_product ( u, w0 )
        e = dot_product ( v, w0 )
    
        !  Check the determinant.
        det = - a * c + b * b
        
        if ( det == 0.0D+00 ) then
            sn = 0.0D+00
            if ( abs ( b ) < abs ( c ) ) then
                tn = e / c
            else
                tn = d / b
            end if
        else
            sn = ( c * d - b * e ) / det
            tn = ( b * d - a * e ) / det
        end if

        !  Now if both nearest points on the lines
        !  also happen to lie inside their line segments,
        !  then we have found the nearest points on the line segments.
        if ( 0.0D+00 <= sn .and. sn <= 1.0D+00 .and. &
          0.0D+00 <= tn .and. tn <= 1.0D+00 ) then
            pn(1:dim_num) = p1(1:dim_num) + sn * ( p2(1:dim_num) - p1(1:dim_num) )
            qn(1:dim_num) = q1(1:dim_num) + tn * ( q2(1:dim_num) - q1(1:dim_num) )
            dist = sqrt ( sum ( ( pn(1:dim_num) - qn(1:dim_num) )**2 ) )
            return
        end if
        !  The nearest point did not occur in the interior.
        !  Therefore it must be achieved at an endpoint.
        !pn, qn, sn, tn 
        !ppn, ssn

        !write(*,*) ' SEGMENT: Endpt contact....'
        ! case qn, p1
        call segment_point_dist_3d ( q1, q2, p1, dist2(1), rn(1, 1:3), un(1) )
        !case qn, p2
        call segment_point_dist_3d ( q1, q2, p2, dist2(2), rn(2, 1:3), un(2) )
        !case pn, q1
        call segment_point_dist_3d ( p1, p2, q1, dist2(3), rn(3, 1:3), un(3) )
        !case pn, q2
        call segment_point_dist_3d ( p1, p2, q2, dist2(4), rn(4, 1:3), un(4) )

        close=1
        do i=2,4
            if (dist2(i) .lt. dist2(close)) close=i
        end do
        
        dist=dist2(close)
        if (close .lt. 3) then
            qn(1:dim_num)=rn(close,1:3)
            tn=un(close)
            if (close .eq. 1) then
                pn(1:dim_num)=p1(1:dim_num)
                sn=0.d0
            else
                pn(1:dim_num)=p2(1:dim_num) 
                sn=1.d0
            end if
        else
            pn(1:dim_num)=rn(close,1:3)
            sn=un(close)
            if (close .eq. 3) then
                qn=q1(1:dim_num)
                tn=0.d0
            else
                qn=q2(1:dim_num) 
                tn=1.d0
            end if
        end if
        
        return
    end

    !*****************************************************************************80
    !
    !! SEGMENT_POINT_DIST_3D: distance ( line segment, point ) in 3D.
    !
        !  Discussion:
        !
        !    A line segment is the finite portion of a line that lies between
        !    two points P1 and P2.
        !
        !    The nearest point will satisfy the condition
        !
        !      PN = (1-T) * P1 + T * P2.
        !
        !    T will always be between 0 and 1.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license. 
        !
        !  Modified:
        !
        !    03 May 2006
        !
        !  Author:
        !
        !    John Burkardt
        !

    subroutine segment_point_dist_3d ( p1, p2, p, dist, pn, t )
        !  Parameters:
        !
        !    Input, real(kind=DBL) P1(3), P2(3), the endpoints of the segment.
        !
        !    Input, real(kind=DBL) P(3), the point whose nearest neighbor on
        !    the line segment is to be determined.
        !
        !    Output, real(kind=DBL) DIST, the distance from the point to the
        !    line segment.
        !
        implicit none
            
        integer ( kind = 4 ), parameter :: dim_num = 3
        real(kind=DBL) bot,t
        real(kind=DBL) dist
        real(kind=DBL),dimension(dim_num) :: p,p1,p2,pn
        !
        !  If the line segment is actually a point, then the answer is easy.
        !
        if ( all ( p1(1:dim_num) == p2(1:dim_num) ) ) then    
            t = 0.0D+00
        else    
            bot = sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 )
            t = sum ( ( p(1:dim_num)  - p1(1:dim_num) ) &
              * ( p2(1:dim_num) - p1(1:dim_num) ) ) / bot
            
            t = max ( t, 0.0D+00 )
            t = min ( t, 1.0D+00 )           
        end if
        pn(1:dim_num) = p1(1:dim_num) + t * ( p2(1:dim_num) - p1(1:dim_num) )
        dist = sqrt ( sum ( ( p(1:dim_num) - pn(1:dim_num) )**2 ) )
        return
    end

end module otter_math