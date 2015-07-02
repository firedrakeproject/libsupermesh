#define BUF_SIZE 24
#include "fdebug.h"

module libsupermesh_tri_intersection_module

  use libsupermesh_fldebug
  use libsupermesh_tri_intersection_2_module, only : tri_type, &
        & intersect_tris_libwm => libsupermesh_intersect_tris_libwm, intersect_tris_dt_old => libsupermesh_intersect_tris_dt_old
  use libsupermesh_tri_intersection_2_module, only : libwm_tri_buf_size => tri_buf_size

  use mpi

  implicit none

  private

  public :: tri_type, intersect_tris, intersect_tris_dt_old, intersect_tris_libwm

  type line_type
    real, dimension(2) :: normal
    real, dimension(2) :: point
  end type line_type

  interface intersect_tris
    module procedure intersect_tris_dt, intersect_tris_dt_public
  end interface intersect_tris

  interface get_lines
    module procedure get_lines_tri
  end interface get_lines

  integer, parameter, public :: tri_buf_size = BUF_SIZE
  public :: libwm_tri_buf_size
  
  real, dimension(2, BUF_SIZE + 2), save :: points_tmp
  integer, save :: n_points_tmp

contains

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
