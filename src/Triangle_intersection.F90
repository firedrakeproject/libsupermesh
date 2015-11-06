#define BUF_SIZE 24
#include "fdebug.h"

module libsupermesh_tri_intersection_module

  use iso_c_binding

  use libsupermesh_fldebug
  use libsupermesh_tri_intersection_2_module, only : tri_type, &
        & intersect_tris_dt_old

  use mpi

  implicit none

  private

  public :: tri_type, intersect_tris, intersect_tris_dt_old

  type line_type
    real, dimension(2) :: normal
    real, dimension(2) :: point
  end type line_type

  interface intersect_tris
    module procedure intersect_tris_dt, intersect_tris_dt_public, &
         & intersect_tris_libwm
  end interface intersect_tris

  interface get_lines
    module procedure get_lines_tri
  end interface get_lines

  integer, parameter, public :: tri_buf_size = BUF_SIZE

  real, dimension(2, BUF_SIZE + 2), save :: points_tmp
  integer, save :: n_points_tmp

  interface cintersector_set_input
    subroutine libsupermesh_cintersector_set_input(nodes_A, nodes_B, ndim, loc) bind(c)
      use iso_c_binding
      implicit none
      integer(kind = c_int), intent(in) :: ndim, loc
      real(kind = c_double), dimension(ndim, loc), intent(in) :: nodes_A, nodes_B
    end subroutine libsupermesh_cintersector_set_input
  end interface cintersector_set_input

  interface cintersector_drive
    subroutine libsupermesh_cintersector_drive() bind(c)
      implicit none
    end subroutine libsupermesh_cintersector_drive
  end interface cintersector_drive

  interface cintersector_query
    subroutine libsupermesh_cintersector_query(nonods, totele) bind(c)
      use iso_c_binding
      implicit none
      integer(kind = c_int), intent(out) :: nonods, totele
    end subroutine libsupermesh_cintersector_query
  end interface cintersector_query

  interface cintersector_get_output
    subroutine libsupermesh_cintersector_get_output(nonods, totele, ndim, loc, nodes, enlist) bind(c)
      use iso_c_binding
      implicit none
      integer(kind = c_int), intent(in) :: nonods, totele, ndim, loc
      real(kind = c_double), dimension(ndim, nonods), intent(out) :: nodes
      integer(kind = c_int), dimension(loc, totele), intent(out) :: enlist
    end subroutine libsupermesh_cintersector_get_output
  end interface cintersector_get_output

  interface cintersector_set_dimension
    subroutine libsupermesh_cintersector_set_dimension(ndim) bind(c)
      use iso_c_binding
      implicit none
      integer(kind = c_int), intent(in) :: ndim
    end subroutine libsupermesh_cintersector_set_dimension
  end interface cintersector_set_dimension

  interface cintersection_finder_reset
    subroutine libsupermesh_cintersection_finder_reset(ntests) bind(c)
      use iso_c_binding
      implicit none
      integer(kind = c_int), intent(out) :: ntests
    end subroutine libsupermesh_cintersection_finder_reset
  end interface cintersection_finder_reset

  public :: cintersector_set_input, cintersector_drive, cintersector_query, &
    & cintersector_get_output, cintersector_set_dimension, &
    & cintersection_finder_reset

contains

  subroutine intersect_tris_libwm(triA, triB, nodesC, ndglnoC, n_trisC)
    real, dimension(2, 3), intent(in) :: triA
    real, dimension(2, 3), intent(in) :: triB
    real, dimension(2, BUF_SIZE), intent(out) :: nodesC
    integer, dimension(3, BUF_SIZE), intent(out) :: ndglnoC
    integer, intent(out) :: n_trisC

    integer :: i, nonods

    call cintersector_set_input(triA, triB, 2, 3)
    call cintersector_drive
    call cintersector_query(nonods, n_trisC)
    call cintersector_get_output(nonods, n_trisC, 2, 3, nodesC, ndglnoC)

  end subroutine intersect_tris_libwm

  subroutine intersect_tris_dt_public(triA, triB, trisC, n_trisC)
    real, dimension(2, 3), intent(in) :: triA
    real, dimension(2, 3), intent(in) :: triB
    real, dimension(2, 3, BUF_SIZE), intent(out) :: trisC
    integer, intent(out) :: n_trisC

    integer :: i
    type(tri_type) :: triA_t, triB_t
    type(tri_type), dimension(BUF_SIZE), save :: trisC_t

    triA_t%V = triA
    triB_t%V = triB
    call intersect_tris_dt(triA_t, triB_t, trisC_t, n_trisC)
    do i = 1, n_trisC
      trisC(:, :, i) = trisC_t(i)%V
    end do

  end subroutine intersect_tris_dt_public

  subroutine intersect_tris_dt(triA, triB, trisC, n_trisC)
    type(tri_type), intent(in) :: triA
    type(tri_type), intent(in) :: triB
    type(tri_type), dimension(BUF_SIZE), intent(out) :: trisC
    integer, intent(out) :: n_trisC

    integer :: i
    real :: tol
    type(line_type), dimension(3) :: lines_b

    real, dimension(2, BUF_SIZE + 2), save :: points
    integer :: n_points

    lines_b = get_lines(triB)

    points(:, :3) = triA%v
    n_points = 3

    n_trisC = 0
    call clip(lines_b(1), points, n_points)
    if(n_points_tmp < 3) return
    points(:, :n_points_tmp) = points_tmp(:, :n_points_tmp)
    n_points = n_points_tmp
    call clip(lines_b(2), points, n_points)
    if(n_points_tmp < 3) return
    points(:, :n_points_tmp) = points_tmp(:, :n_points_tmp)
    n_points = n_points_tmp
    call clip(lines_b(3), points, n_points)
    if(n_points_tmp < 3) return

    tol = 10.0 * min(spacing(triangle_area(triA)), spacing(triangle_area(triB)))
    do i = 1, n_points_tmp - 2
      n_trisC = n_trisC + 1
      trisC(n_trisC)%V(:, 1) = points_tmp(:, 1)
      trisC(n_trisC)%V(:, 2) = points_tmp(:, i + 1)
      trisC(n_trisC)%V(:, 3) = points_tmp(:, i + 2)
      if(triangle_area(trisC(n_trisC)) < tol) then
        n_trisC = n_trisC - 1
      end if
    end do

  end subroutine intersect_tris_dt

  ! Sutherland-Hodgman clipping algorithm
  subroutine clip(line, points, n_points)
    type(line_type), intent(in) :: line
    real, dimension(2, BUF_SIZE + 2), intent(in) :: points
    integer, intent(in) :: n_points

    integer :: i
    real :: d1, d2, f
    real, dimension(2) :: p1, p2
    real, dimension(BUF_SIZE + 2), save :: d

    do i = 1, n_points
      d(i) = dot_product(line%normal, points(:, i) - line%point)
    end do

    n_points_tmp = 0
    do i = 1, n_points
      p1 = points(:, i)
      d1 = d(i)
      if(i == n_points) then
        p2 = points(:, 1)
        d2 = d(1)
      else
        p2 = points(:, i + 1)
        d2 = d(i + 1)
      end if

      if(d1 <= 0.0) then
        if(d2 <= 0.0) then
          ! No clip
          n_points_tmp = n_points_tmp + 1
          points_tmp(:, n_points_tmp) = p1
        else
          ! New point
          n_points_tmp = n_points_tmp + 1
          points_tmp(:, n_points_tmp) = p1
          n_points_tmp = n_points_tmp + 1
          f = max(min((abs(d1) / (abs(d1) + abs(d2))), 1.0), 0.0)
          points_tmp(:, n_points_tmp) = p1 + f * (p2 - p1)
        end if
      else if(d2 <= 0.0) then
        ! Move point
        n_points_tmp = n_points_tmp + 1
        f = max(min((abs(d1) / (abs(d1) + abs(d2))), 1.0), 0.0)
        points_tmp(:, n_points_tmp) = p1 + f * (p2 - p1)
      !else
        ! Full clip
      end if
    end do

  end subroutine clip

  pure function get_lines_tri(tri) result(lines)
    type(tri_type), intent(in) :: tri
    type(line_type), dimension(3) :: lines

    real :: det

    ! Note that the normals are not normalised

    lines(1)%normal(1) = -(tri%V(2, 2) - tri%V(2, 1))
    lines(1)%normal(2) =  (tri%V(1, 2) - tri%V(1, 1))
    lines(1)%point = tri%V(:, 1)

    lines(2)%normal(1) = -(tri%V(2, 3) - tri%V(2, 2))
    lines(2)%normal(2) =  (tri%V(1, 3) - tri%V(1, 2))
    lines(2)%point = tri%V(:, 2)

    lines(3)%normal(1) = -(tri%V(2, 1) - tri%V(2, 3))
    lines(3)%normal(2) =  (tri%V(1, 1) - tri%V(1, 3))
    lines(3)%point = tri%V(:, 3)

    ! Transform the normals so that their direction is consistent (step b).
    det = dot_product(tri%V(:, 2) - tri%V(:, 1), lines(3)%normal)
    if(det > 0.0) then
      lines(1)%normal = -lines(1)%normal
      lines(2)%normal = -lines(2)%normal
      lines(3)%normal = -lines(3)%normal
    end if

  end function get_lines_tri

  pure function triangle_area(tri) result(area)
    type(tri_type), intent(in) :: tri

    real :: area
    real, dimension(2) :: u, v

    u = tri%V(:, 3) - tri%V(:, 1)
    v = tri%V(:, 2) - tri%V(:, 1)

    area = 0.5 * abs(u(2) * v(1) - u(1) * v(2))

  end function triangle_area

end module libsupermesh_tri_intersection_module
