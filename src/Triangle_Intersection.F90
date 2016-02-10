#include "fdebug.h"

#define BUF_SIZE 22

module libsupermesh_tri_intersection

  implicit none

  private

  public :: tri_type, line_type, max_n_trisC, intersect_tris, intersect_polys, &
    & get_lines, triangle_area

  type tri_type
    real, dimension(2, 3) :: v
  end type tri_type

  type line_type
    real, dimension(2) :: normal
    real, dimension(2) :: point
  end type line_type
  
  interface max_n_trisC
    module procedure max_n_trisC_tri, max_n_trisC_poly
  end interface max_n_trisC

  interface intersect_tris
    module procedure intersect_tris_real, intersect_tris_tri
  end interface intersect_tris
  
  interface intersect_polys
    module procedure intersect_tri_lines, intersect_polys_real, &
      & intersect_polys_lines
  end interface intersect_polys

  interface get_lines
    module procedure get_lines_tri, get_lines_poly
  end interface get_lines
  
  interface triangle_area
    module procedure triangle_area_real, triangle_area_tri
  end interface triangle_area

  integer, parameter, public :: tri_buf_size = BUF_SIZE

  real, dimension(2, BUF_SIZE + 2), save :: points_tmp
  integer, save :: n_points_tmp

contains

  pure function max_n_trisC_tri(n_linesB) result(max_n_trisC)
    integer, intent(in) :: n_linesB
    
    integer :: max_n_trisC
    
    max_n_trisC = 3 * (2 ** n_linesB) - 2
    
  end function max_n_trisC_tri

  pure function max_n_trisC_poly(n_linesA, n_linesB) result(max_n_trisC)
    integer, intent(in) :: n_linesA
    integer, intent(in) :: n_linesB
    
    integer :: max_n_trisC
    
    max_n_trisC = n_linesA * (2 ** n_linesB) - 2
    
  end function max_n_trisC_poly

  subroutine intersect_tris_real(triA, triB, trisC, n_trisC)
    real, dimension(2, 3), intent(in) :: triA
    real, dimension(2, 3), intent(in) :: triB
    real, dimension(2, 3, BUF_SIZE), intent(out) :: trisC
    integer, intent(out) :: n_trisC

    integer :: i
    type(tri_type) :: triA_t, triB_t
    type(tri_type), dimension(BUF_SIZE), save :: trisC_t

    triA_t%v = triA
    triB_t%v = triB
    call intersect_tris(triA_t, triB_t, trisC_t, n_trisC)
    do i = 1, n_trisC
      trisC(:, :, i) = trisC_t(i)%v
    end do

  end subroutine intersect_tris_real

  subroutine intersect_tris_tri(triA, triB, trisC, n_trisC)
    type(tri_type), intent(in) :: triA
    type(tri_type), intent(in) :: triB
    type(tri_type), dimension(BUF_SIZE), intent(out) :: trisC
    integer, intent(out) :: n_trisC

    integer :: i
    real :: tol
    type(line_type), dimension(3) :: linesB

    real, dimension(2, BUF_SIZE + 2), save :: points
    integer :: n_points

    linesB = get_lines(triB)

    points(:, :3) = triA%v
    n_points = 3

    n_trisC = 0
    call clip_buf(linesB(1), points, n_points)
    if(n_points_tmp < 3) return
    points(:, :n_points_tmp) = points_tmp(:, :n_points_tmp)
    n_points = n_points_tmp
    call clip_buf(linesB(2), points, n_points)
    if(n_points_tmp < 3) return
    points(:, :n_points_tmp) = points_tmp(:, :n_points_tmp)
    n_points = n_points_tmp
    call clip_buf(linesB(3), points, n_points)
    if(n_points_tmp < 3) return

    tol = 10.0 * min(spacing(triangle_area(triA)), spacing(triangle_area(triB)))
    do i = 1, n_points_tmp - 2
      n_trisC = n_trisC + 1
      trisC(n_trisC)%v(:, 1) = points_tmp(:, 1)
      trisC(n_trisC)%v(:, 2) = points_tmp(:, i + 1)
      trisC(n_trisC)%v(:, 3) = points_tmp(:, i + 2)
      if(triangle_area(trisC(n_trisC)) < tol) then
        n_trisC = n_trisC - 1
      end if
    end do

  end subroutine intersect_tris_tri

  subroutine intersect_tri_lines(triA, linesB, trisC, n_trisC, areaB, work)
    type(tri_type), intent(in) :: triA
    type(line_type), dimension(:), intent(in) :: linesB
    type(tri_type), dimension(:), intent(out) :: trisC
    integer, intent(out) :: n_trisC
    real, optional, intent(in) :: areaB
    real, dimension(:, :, :), target, optional, intent(out) :: work

    call intersect_polys(triA%v, linesB, trisC, n_trisC, areaA = triangle_area(triA), areaB = areaB, work = work)

  end subroutine intersect_tri_lines
  
  subroutine intersect_polys_real(polyA, polyB, trisC, n_trisC, areaA, areaB, work)
    ! 2 x loc_a
    real, dimension(:, :), intent(in) :: polyA
    ! 2 x loc_b
    real, dimension(:, :), intent(in) :: polyB
    type(tri_type), dimension(:), intent(out) :: trisC
    integer, intent(out) :: n_trisC
    real, optional, intent(in) :: areaA
    real, optional, intent(in) :: areaB
    real, dimension(:, :, :), target, optional, intent(out) :: work
    
    call intersect_polys(polyA, get_lines(polyB), trisC, n_trisC, areaA = areaA, areaB = areaB, work = work)
  
  end subroutine intersect_polys_real

  subroutine intersect_polys_lines(polyA, linesB, trisC, n_trisC, areaA, areaB, work)
    ! 2 x loc_a
    real, dimension(:, :), intent(in) :: polyA
    ! loc_b
    type(line_type), dimension(:), intent(in) :: linesB
    type(tri_type), dimension(:), intent(out) :: trisC
    integer, intent(out) :: n_trisC
    real, optional, intent(in) :: areaA
    real, optional, intent(in) :: areaB
    real, dimension(:, :, :), target, optional, intent(out) :: work

    integer :: i
    real :: tol

    real, dimension(:, :), pointer :: points, points_new
    integer :: n_points, n_points_new
    
    if(present(work)) then
      points => work(:, :, 1)
      points_new => work(:, :, 2)
    else
      allocate(points(2, max_n_trisC(size(polyA, 2), size(linesB)) + 2))
      allocate(points_new(2, size(points)))
    end if

    points(:, :size(polyA, 2)) = polyA
    n_points = size(polyA, 2)

    n_trisC = 0
    do i = 1, size(linesB)
      call clip(linesB(i), points, n_points, points_new, n_points_new)
      if(n_points_new < 3) goto 42
      points(:, :n_points_new) = points_new(:, :n_points_new)
      n_points = n_points_new
    end do

    if(present(areaB)) then
      if(present(areaA)) then
        tol = 10.0 * min(spacing(areaA), spacing(areaB))
      else
        tol = 10.0 * spacing(areaB)
      end if
    else if(present(areaA)) then
      tol = 10.0 * spacing(areaA)
    else
      tol = 10.0 * epsilon(0.0)
    end if
    do i = 1, n_points_new - 2
      n_trisC = n_trisC + 1
      trisC(n_trisC)%v(:, 1) = points_new(:, 1)
      trisC(n_trisC)%v(:, 2) = points_new(:, i + 1)
      trisC(n_trisC)%v(:, 3) = points_new(:, i + 2)
      if(triangle_area(trisC(n_trisC)) < tol) then
        n_trisC = n_trisC - 1
      end if
    end do

42  if(.not. present(work)) deallocate(points, points_new)

  end subroutine intersect_polys_lines

  ! Sutherland-Hodgman clipping algorithm. See:
  !   Reentrant polygon clipping, I. E. Sutherland and G. W. Hodgman,
  !   Communications of the ACM vol. 17, 1974, pp. 32--42.
  subroutine clip_buf(line, points, n_points)
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

  end subroutine clip_buf

  ! Sutherland-Hodgman clipping algorithm. See:
  !   Reentrant polygon clipping, I. E. Sutherland and G. W. Hodgman,
  !   Communications of the ACM vol. 17, 1974, pp. 32--42.
  subroutine clip(line, points, n_points, points_new, n_points_new)
    type(line_type), intent(in) :: line
    real, dimension(:, :), intent(in) :: points
    integer, intent(in) :: n_points
    real, dimension(:, :), intent(out) :: points_new
    integer, intent(out) :: n_points_new

    integer :: i
    real :: d0, d1, d2, f
    real, dimension(2) :: p1, p2

    d0 = dot_product(line%normal, points(:, 1) - line%point)
    d2 = d0
    n_points_new = 0
    do i = 1, n_points
      p1 = points(:, i)
      d1 = d2
      if(i == n_points) then
        p2 = points(:, 1)
        d2 = d0
      else
        p2 = points(:, i + 1)
        d2 = dot_product(line%normal, points(:, i + 1) - line%point)
      end if

      if(d1 <= 0.0) then
        if(d2 <= 0.0) then
          ! No clip
          n_points_new = n_points_new + 1
          points_new(:, n_points_new) = p1
        else
          ! New point
          n_points_new = n_points_new + 1
          points_new(:, n_points_new) = p1
          n_points_new = n_points_new + 1
          f = max(min((abs(d1) / (abs(d1) + abs(d2))), 1.0), 0.0)
          points_new(:, n_points_new) = p1 + f * (p2 - p1)
        end if
      else if(d2 <= 0.0) then
        ! Move point
        n_points_new = n_points_new + 1
        f = max(min((abs(d1) / (abs(d1) + abs(d2))), 1.0), 0.0)
        points_new(:, n_points_new) = p1 + f * (p2 - p1)
      !else
        ! Full clip
      end if
    end do

  end subroutine clip

  pure function get_lines_tri(tri) result(lines)
    type(tri_type), intent(in) :: tri
    
    type(line_type), dimension(3) :: lines

    ! Note that the normals are not normalised

    lines(1)%normal(1) = -(tri%v(2, 2) - tri%v(2, 1))
    lines(1)%normal(2) =  (tri%v(1, 2) - tri%v(1, 1))
    lines(1)%point = tri%v(:, 1)

    lines(2)%normal(1) = -(tri%v(2, 3) - tri%v(2, 2))
    lines(2)%normal(2) =  (tri%v(1, 3) - tri%v(1, 2))
    lines(2)%point = tri%v(:, 2)

    lines(3)%normal(1) = -(tri%v(2, 1) - tri%v(2, 3))
    lines(3)%normal(2) =  (tri%v(1, 1) - tri%v(1, 3))
    lines(3)%point = tri%v(:, 3)

    if(dot_product(tri%v(:, 2) - tri%v(:, 1), lines(3)%normal) > 0.0) then
      lines(1)%normal = -lines(1)%normal
      lines(2)%normal = -lines(2)%normal
      lines(3)%normal = -lines(3)%normal
    end if

  end function get_lines_tri

  pure function get_lines_poly(poly) result(lines)
    ! 2 x loc
    real, dimension(:, :), intent(in) :: poly
    
    type(line_type), dimension(size(poly, 2)) :: lines

    integer :: i, loc
    
    loc = size(poly, 2)
  
    ! Assumes clockwise or anti-clockwise ordering
    ! Note that the normals are not normalised

    do i = 1, loc - 1
      lines(i)%normal(1) =  (poly(2, i + 1) - poly(2, i))
      lines(i)%normal(2) = -(poly(1, i + 1) - poly(1, i))
      lines(i)%point = poly(:, i)
    end do
    lines(loc)%normal(1) =  (poly(2, 1) - poly(2, loc))
    lines(loc)%normal(2) = -(poly(1, 1) - poly(1, loc))
    lines(loc)%point = poly(:, loc)
    
    if(dot_product(poly(:, 1) - lines(2)%point, lines(2)%normal) > 0.0) then
      do i = 1, loc
        lines(i)%normal = -lines(i)%normal
      end do
    end if
 
  end function get_lines_poly

  pure function triangle_area_real(tri) result(area)
    real, dimension(2, 3), intent(in) :: tri

    real :: area
    real, dimension(2) :: u, v

    u = tri(:, 3) - tri(:, 1)
    v = tri(:, 2) - tri(:, 1)

    area = 0.5 * abs(u(2) * v(1) - u(1) * v(2))

  end function triangle_area_real

  pure function triangle_area_tri(tri) result(area)
    type(tri_type), intent(in) :: tri

    real :: area
    real, dimension(2) :: u, v

    u = tri%v(:, 3) - tri%v(:, 1)
    v = tri%v(:, 2) - tri%v(:, 1)

    area = 0.5 * abs(u(2) * v(1) - u(1) * v(2))

  end function triangle_area_tri
  
end module libsupermesh_tri_intersection