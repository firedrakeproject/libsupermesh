#include "fdebug.h"

subroutine test_tri_intersector()

  use libsupermesh_intersection_finder, only : intersections, deallocate, &
    & intersection_finder
  use libsupermesh_read_triangle, only : read_ele, read_node
  use libsupermesh_supermesh, only : intersect_elements, &
    & intersect_simplices, triangle_area
  use libsupermesh_tri_intersection, only : tri_type, tri_buf_size, &
    & intersect_tris, intersect_polys, get_lines
  use libsupermesh_unittest_tools, only : report_test, operator(.fne.)
  
  implicit none

  integer :: ele_a, ele_b, ele_c, i, n_trisC
  integer, dimension(:, :), allocatable :: enlist_a, enlist_b
  real :: area_c
  real, dimension(2, 3) :: triA_real
  real, dimension(2, 3, tri_buf_size) :: trisC_real
  real, dimension(2, tri_buf_size + 2, 2) :: work
  real, dimension(:, :), allocatable :: positions_a, positions_b  
  type(intersections), dimension(:), allocatable :: map_ab
  type(tri_type) :: triA, triB
  type(tri_type), dimension(tri_buf_size) :: trisC

  integer, parameter :: dim = 2
  
  call read_node("data/triangle_0_01.node", dim, positions_a)
  call read_ele("data/triangle_0_01.ele", dim, enlist_a)
  call read_node("data/square_0_01.node", dim, positions_b)
  call read_ele("data/square_0_01.ele", dim, enlist_b)
  
  allocate(map_ab(size(enlist_a, 2)))
  call intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)

  area_c = 0.0
  do ele_a = 1, size(enlist_a, 2)
    triA%v = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      triB%v = positions_b(:, enlist_b(:, ele_b))
      call intersect_tris(triA, triB, trisC, n_trisC)
      do ele_c = 1, n_trisC
        area_c = area_c + triangle_area(trisC(ele_c)%v)
      end do
    end do    
  end do
  call report_test("[intersect_tris]", area_c .fne. 0.5, .false., "Incorrect intersection area")

  area_c = 0.0
  do ele_a = 1, size(enlist_a, 2)
    triA_real = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_tris(triA_real, positions_b(:, enlist_b(:, ele_b)), trisC_real, n_trisC)
      do ele_c = 1, n_trisC
        area_c = area_c + triangle_area(trisC_real(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_tris]", area_c .fne. 0.5, .false., "Incorrect intersection area")

  area_c = 0.0
  do ele_a = 1, size(enlist_a, 2)
    triA%v = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      triB%v = positions_b(:, enlist_b(:, ele_b))
      call intersect_polys(triA, get_lines(triB), trisC, n_trisC, work = work)
      do ele_c = 1, n_trisC
        area_c = area_c + triangle_area(trisC(ele_c)%v)
      end do
    end do    
  end do
  call report_test("[intersect_polys]", area_c .fne. 0.5, .false., "Incorrect intersection area")

  area_c = 0.0
  do ele_a = 1, size(enlist_a, 2)
    triA_real = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_polys(triA_real, get_lines(positions_b(:, enlist_b(:, ele_b))), trisC, n_trisC, work = work)
      do ele_c = 1, n_trisC
        area_c = area_c + triangle_area(trisC(ele_c)%v)
      end do
    end do    
  end do
  call report_test("[intersect_polys]", area_c .fne. 0.5, .false., "Incorrect intersection area")

  area_c = 0.0
  do ele_a = 1, size(enlist_a, 2)
    triA_real = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_polys(triA_real, positions_b(:, enlist_b(:, ele_b)), trisC, n_trisC, work = work)
      do ele_c = 1, n_trisC
        area_c = area_c + triangle_area(trisC(ele_c)%v)
      end do
    end do    
  end do
  call report_test("[intersect_polys]", area_c .fne. 0.5, .false., "Incorrect intersection area")
  
  area_c = 0.0
  do ele_a = 1, size(enlist_a, 2)
    triA_real = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_simplices(triA_real, positions_b(:, enlist_b(:, ele_b)), trisC_real, n_trisC)
      do ele_c = 1, n_trisC
        area_c = area_c + triangle_area(trisC_real(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_simplices]", area_c .fne. 0.5, .false., "Incorrect intersection area")
  
  area_c = 0.0
  do ele_a = 1, size(enlist_a, 2)
    triA_real = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_elements(triA_real, positions_b(:, enlist_b(:, ele_b)), trisC_real, n_trisC)
      do ele_c = 1, n_trisC
        area_c = area_c + triangle_area(trisC_real(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_elements]", area_c .fne. 0.5, .false., "Incorrect intersection area")
  
  call deallocate(map_ab)
  deallocate(map_ab, positions_a, enlist_a, positions_b, enlist_b)

end subroutine test_tri_intersector
