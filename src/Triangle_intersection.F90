#define BUF_SIZE 8
#include "fdebug.h"

!
! Original Triangle (A)
!
!         E(Ex, Ey)
!         /\
!        /  \
!       /    \
!      /      \
!     /--------\
!  N(Nx,Ny)     D(Dx,Dy)
!
! <==================================>
!
! Clipping triangle (B)
!
!       Q(Qx, Qy)
!       /\
!      /  \
!     /    \
!    /      \
!   /--------\
! P(Px,Py)    R(Rx, Ry)
!
! <==================================>
!
! Intersecting Triangles
!
!      Q(Qx,Qy)  E(Ex, Ey)
!            |   /\
!  <---------|  /  \
!  ~         | /    \
!  n         |/      \
!            *        \
!           /|         \
!          / |          \
!         /  |           \
!        /   |            \
!       /    |             \
!      /     |              \
!     /      |I(Ix,Iy)       \
!    /-------*----------------\
!   N(Nx,Ny) |            D(Dx,Dy)
!            |
!            |
!O(Ox,Oy)    P(Px,Py)
!
! Steps:
! Sutherland–Hodgman Algorithm
! a. Find the normals of all sides of Triangle B.
! b. Transform the normals so that their direction is consistent.
!                                                    ->   ~
! c. Decide if we want to keep a node of Triangle A (PN * n >0).
!                                ->    ->   ->       ->         ->   ->
! d. Find the co-ordinates of I (NI = -ON + OI = l * ND = l * (-ON + OD)


module libsupermesh_tri_intersection_2_module

  use libsupermesh_fldebug

  implicit none

  type tri_type
    real, dimension(2, 3) :: V ! vertices of the tri
  end type tri_type

  type line_type
    real, dimension(2) :: normal
    real :: c = 0
  end type line_type

  type(tri_type), dimension(BUF_SIZE), save :: tri_array_tmp
  integer :: tri_cnt_tmp = 0

  public :: tri_type, line_type, get_lines, libsupermesh_intersect_tris_dt_old

  interface get_lines
    module procedure get_lines_tri
  end interface

  integer, parameter, public :: tri_buf_size = BUF_SIZE

contains

  subroutine libsupermesh_intersect_tris_dt_old(triA, triB, trisC, n_trisC)
    type(tri_type), intent(in) :: triA
    type(tri_type), intent(in) :: triB
    type(tri_type), dimension(BUF_SIZE), intent(out) :: trisC
    integer, intent(out) :: n_trisC

    integer :: i, j
    type(line_type), dimension(3) :: linesB

    linesB = get_lines(triB)
    n_trisC = 1
    trisC(1) = triA
    do i=1,size(linesB)
      ! Clip the trisC against the i'th plane
      tri_cnt_tmp = 0

      do j=1,n_trisC
        call clip(linesB(i), trisC(j))
      end do

      if (i /= size(linesB)) then
        n_trisC = tri_cnt_tmp
        trisC(1:n_trisC) = tri_array_tmp(1:n_trisC)
      else
        ! Copy the result
        n_trisC = 0
        do j=1,tri_cnt_tmp
            n_trisC = n_trisC + 1
            trisC(n_trisC) = tri_array_tmp(j)
        end do
      end if
    end do
  end subroutine libsupermesh_intersect_tris_dt_old
  
  subroutine clip(line, tri)
  ! Clip tri against the plane
  ! and append any output to tri_array_tmp.
    type(line_type), intent(in) :: line
    type(tri_type), intent(in) :: tri

    real, dimension(3) :: dists
    integer :: neg_cnt, pos_cnt
    integer, dimension(3) :: neg_idx, pos_idx
    integer :: i

    real :: invdiff, w0, w1
    type(tri_type) :: tri_tmp

    ! Negative == inside
    ! Positive == outside

    neg_cnt = 0
    pos_cnt = 0
    dists = distances_to_line(line, tri)

    do i=1,3
      if (dists(i) <= 0.0) then
        neg_cnt = neg_cnt + 1
        neg_idx(neg_cnt) = i
      else
        pos_cnt = pos_cnt + 1
        pos_idx(pos_cnt) = i
      end if
    end do

    if (pos_cnt == 0) then
      ! tri is completely on positive side of line, full clip
      return
    end if

    if (neg_cnt == 0) then
      ! tri is completely on negative side of line, no clip
      tri_cnt_tmp = tri_cnt_tmp + 1
      tri_array_tmp(tri_cnt_tmp) = tri
      return
    end if

    ! The tri is split by the line, so we have more work to do.

    select case(pos_cnt)
    case(1)
      select case(neg_cnt)
      case(2)
        tri_cnt_tmp = tri_cnt_tmp + 1
        tri_array_tmp(tri_cnt_tmp) = tri

        do i=1,neg_cnt
          invdiff = 1.0 / ( dists(pos_idx(1)) - dists(neg_idx(i)) )
          w0 = -dists(neg_idx(i)) * invdiff
          w1 =  dists(pos_idx(1)) * invdiff
          tri_array_tmp(tri_cnt_tmp)%V(:, neg_idx(i)) = &
               w0 * tri%V(:, pos_idx(1)) + w1 * tri%V(:, neg_idx(i))
        end do
      case default
        FLAbort("Error. Found more than three points.")
      end select
    case(2)
      select case(neg_cnt)
      case(1)
        ! We need to return two triangles back
        do i=1,pos_cnt
          invdiff = 1.0 / ( dists(pos_idx(i)) - dists(neg_idx(1)) )
          w0 = -dists(neg_idx(1)) * invdiff
          w1 =  dists(pos_idx(i)) * invdiff
          tri_tmp%V(:,i) = &
               w0 * tri%V(:, pos_idx(i)) + w1 * tri%V(:, neg_idx(1))
        end do

        tri_cnt_tmp = tri_cnt_tmp + 1
        tri_array_tmp(tri_cnt_tmp) = tri
        tri_array_tmp(tri_cnt_tmp)%V(:, pos_idx(1)) = tri_tmp%V(:, 1)
        tri_array_tmp(tri_cnt_tmp)%V(:, pos_idx(2)) = tri_tmp%V(:, 2)
        tri_array_tmp(tri_cnt_tmp)%V(:, neg_idx(1)) = tri%V(:, pos_idx(1))

        tri_cnt_tmp = tri_cnt_tmp + 1
        tri_array_tmp(tri_cnt_tmp) = tri
        tri_array_tmp(tri_cnt_tmp)%V(:, pos_idx(1)) = tri%V(:, pos_idx(1))
        tri_array_tmp(tri_cnt_tmp)%V(:, pos_idx(2)) = tri_tmp%V(:, 2)
        tri_array_tmp(tri_cnt_tmp)%V(:, neg_idx(1)) = tri%V(:, pos_idx(2))
      case default
        FLAbort("Error. Found more than three points.")
      end select
    end select
  end subroutine clip

! Find the normals of all sides of triangle B (steps a and b).
!  pure function get_lines_tri(tri) result(lines)
  function get_lines_tri(tri) result(lines)
    type(tri_type), intent(in) :: tri
    type(line_type), dimension(3) :: lines

    real, dimension(2) :: edge10, edge20, edge30
    real :: det
    integer :: i
    
    ! Find the normals of all sides of triangle B (step a).
    ! P = 1, Q = 2, R = 3
    edge10 = tri%V(:, 2) - tri%V(:, 1); ! PQ = ( Qx - Px , Qy - Py )
    edge20 = tri%V(:, 3) - tri%V(:, 2); ! QR = ( Rx - Qx , Ry - Qy )
    edge30 = tri%V(:, 1) - tri%V(:, 3); ! RP = ( Px - Rx , Py - Ry )

    lines(1)%normal = unit_cross(edge10) ! PQ normal = ( -PQy , PQx )
    lines(2)%normal = unit_cross(edge20) ! QR normal = ( -QRy , QRx )
    lines(3)%normal = unit_cross(edge30) ! PR normal = ( -PRy , PRx )

    ! Transform the normals so that their direction is consistent (step b).
    det = dot_product(edge10, lines(3)%normal)
    if (det < 0) then
      do i=1,3
        lines(i)%normal = -lines(i)%normal
      end do
    end if

    ! And calibrate what is the zero of this line by dotting with
    ! a point we know to be on it
    do i=1,3
      lines(i)%c = dot_product(tri%V(:, i), lines(i)%normal)
    end do

  end function get_lines_tri

! Return a cross product (in 2D).  
  pure function unit_cross(vecA) result(cross)
    real, dimension(2), intent(in) :: vecA
    real, dimension(2) :: cross
    cross(1) = (-1) * vecA(2)
    cross(2) = vecA(1)

  end function unit_cross

  pure function distances_to_line(line, tri) result(dists)
    type(line_type), intent(in) :: line
    type(tri_type), intent(in) :: tri
    real, dimension(3) :: dists
    integer :: i

    forall(i=1:3)
      dists(i) = dot_product(line%normal, tri%V(:, i)) - line%c
    end forall

  end function distances_to_line

end module libsupermesh_tri_intersection_2_module
