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


module libsupermesh_tri_intersection_module

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

  interface libsupermesh_intersect_tris
    module procedure libsupermesh_intersect_tris_dt, libsupermesh_intersect_tris_dt_public, &
     &  libsupermesh_intersect_tris_libwm
  end interface

  interface get_lines
    module procedure get_lines_tri
  end interface

  integer, parameter, public :: tri_buf_size = BUF_SIZE

contains

  subroutine libsupermesh_intersect_tris_libwm(triA, triB, nodesC, ndglnoC, n_trisC)
    real, dimension(2, 3), intent(in) :: triA
    real, dimension(2, 3), intent(in) :: triB
    real, dimension(2, BUF_SIZE), intent(out) :: nodesC
    integer, dimension(3, BUF_SIZE), intent(out) :: ndglnoC
    integer, intent(out) :: n_trisC

    integer :: i, nonods
    real, dimension(2 * BUF_SIZE), save :: lnodesC = -huge(0.0)
    
    call libsupermesh_cintersector_set_input(triA, triB, 2, 3)
    call libsupermesh_cintersector_drive
    call libsupermesh_cintersector_query(nonods, n_trisC)
    call libsupermesh_cintersector_get_output(nonods, n_trisC, 2, 3, lnodesC, ndglnoC)

    do i = 1, nonods
      nodesC(1, i) = lnodesC(i)
      nodesC(2, i) = lnodesC(i + nonods)
    end do

  end subroutine libsupermesh_intersect_tris_libwm

  subroutine libsupermesh_intersect_tris_dt_public(triA, triB, trisC, n_trisC)
    real, dimension(2, 3), intent(in) :: triA
    real, dimension(2, 3), intent(in) :: triB
    real, dimension(2, 3, BUF_SIZE), intent(out) :: trisC
    integer, intent(out) :: n_trisC

    integer :: i
    type(tri_type) :: triA_t, triB_t
    type(tri_type), dimension(BUF_SIZE) :: trisC_t

    triA_t%v = triA
    triB_t%v = triB
    call libsupermesh_intersect_tris_dt(triA_t, triB_t, trisC_t, n_trisC)
    do i = 1, n_trisC
      trisC(:, :, i) = trisC_t(i)%v
    end do

  end subroutine libsupermesh_intersect_tris_dt_public

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
  
  subroutine libsupermesh_intersect_tris_dt(triA, triB, trisC, n_trisC)
    type(tri_type), intent(in) :: triA
    type(tri_type), intent(in) :: triB
    type(tri_type), dimension(BUF_SIZE), intent(out) :: trisC
    integer, intent(out) :: n_trisC

    integer :: i, j, dir, length
    real :: area
    real, dimension(2, 8) :: p1 = 0, p2 = 0, tmp = 0
    type(tri_type) :: temp_tri

    n_trisC = 0
    ! triB == sub       triA == clip
    length = 3
    dir = poly_winding(triA)

    call poly_edge_clip(triB%v, triA%v(:,3), triA%v(:,1), dir, length, p2)
    do i=1,2
      tmp = p2
      p2 = p1
      p1 = tmp
      if ( length .le. 0 ) then 
        length = 1
      end if
      call poly_edge_clip(p1, triA%v(:,i), triA%v(:,i+1), dir, length, p2)
    end do

    do i=2, length-1, 1
      temp_tri%v(:,1) = p2(:, 1)
      temp_tri%v(:,2) = p2(:, i)
      temp_tri%v(:,3) = p2(:, i+1)
      
!      area = triangle_area(temp_tri%v)
!      if ( area .gt. epsilon(0.0) ) then
        n_trisC = n_trisC + 1
        trisC(n_trisC) = temp_tri
!      end if
    end do

contains

  pure function triangle_area(tri) result(area)
    real, dimension(2, 3), intent(in) :: tri

    real :: area
    real, dimension(2) :: u, v

    u = tri(:, 3) - tri(:, 1)
    v = tri(:, 2) - tri(:, 1)

    area = 0.5 * abs(u(2) * v(1) - u(1) * v(2))

  end function triangle_area

  end subroutine libsupermesh_intersect_tris_dt
  
  function poly_winding(tri) result(direction)
    type(tri_type), intent(in) :: tri
    integer :: direction

    direction = left_of(tri%v(:,1), tri%v(:,2), tri%v(:,3))
    
    if ( direction .eq. 0 ) then
      direction = -1
    end if
  end function poly_winding
  
  function left_of(v1, v2, v3) result(direction)
    real, dimension(2), intent(in) :: v1, v2, v3
    integer :: direction

    real, dimension(2) :: tmp1, tmp2
    real :: cross

    tmp1 = vsub(v2, v1)
    tmp2 = vsub(v3, v2)
    cross = tmp1(1) * tmp2(2) - tmp1(2) * tmp2(1)

    if ( cross .lt. -epsilon(0.0) ) then ! HACK ToDO TODO FIX THIS
      direction = -1
    else if ( cross .gt. epsilon(0.0) ) then
      direction = 1
    else
      direction = 0
    end if
  end function left_of
  
  subroutine poly_edge_clip(sub, x0, x1, left, length, p2)
    real, dimension(:, :), intent(in) :: sub
    real, dimension(2), intent(in) :: x0, x1
    integer, intent(in) :: left
    integer, intent(inout) :: length
    real, dimension(:, :), intent(out) :: p2

    real, dimension(2) :: v0, v1, tmp
    real :: cross
    integer :: side0, side1, i, stat, counter

    counter = 1
    v0 = sub(:, length)
    side0 = left_of(x0, x1, v0)

    if ( side0 .ne. -left ) then
      p2(:, counter) = v0
      counter = counter + 1
    end if
    
    do i=1,length
      v1 = sub(:, i)
      side1 = left_of(x0, x1, v1)

      if ( (side0 + side1 .eq. 0) .and. (side0 .ne. 0) ) then
        ! last point and current straddle the edge
        call line_sect(x0, x1, v0, v1, tmp, stat)
        if ( stat .eq. 1 ) then

          p2(:, counter) = tmp
          counter = counter + 1
        end if
      end if
      
      if (i .eq. length ) exit
      if ( side1 .ne. -left ) then

        p2(:, counter) = v1
        counter = counter + 1
      end if
      v0 = v1
      side0 = side1
    end do
    length = counter - 1

  end subroutine poly_edge_clip
  
  subroutine line_sect(x0, x1, y0, y1, res, stat)
    real, dimension(2), intent(in) :: x0, x1, y0, y1
    real, dimension(2), intent(out) :: res
    integer, intent(out) :: stat
    
    real, dimension(2) :: dx, dy, d
    real :: dyx

    stat = 0
    dx = vsub(x1, x0)
    dy = vsub(y1, y0)
    d  = vsub(x0, y0)
!   x0 + a dx = y0 + b dy ->
!   x0 X dx = y0 X dx + b dy X dx ->
!   b = (x0 - y0) X dx / (dy X dx)
    dyx = dy(1) * dx(2) - dy(2) * dx(1)

    if ( abs(dyx) < epsilon(0.0) ) then
      return
    end if

    dyx = ( d(1) * dx(2) - d(2) * dx(1)) / dyx

    if ( ( dyx .le. epsilon(0.0) ) .OR. ( dyx .ge. 1.0 ) ) then
      return
    end if

    res(1) = y0(1) + dyx * dy(1)
    res(2) = y0(2) + dyx * dy(2)

    stat = 1
  
  end subroutine line_sect
  
  function vsub(v1, v2) result(v3)
    real, dimension(2), intent(in) :: v1, v2
    real, dimension(2) :: v3

    v3(1) = v1(1) - v2(1)
    v3(2) = v1(2) - v2(2)
  end function vsub

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

end module libsupermesh_tri_intersection_module
