#include "libsupermesh_debug.h"

subroutine test_intersection_finder_3d()

  use libsupermesh_intersection_finder, only : intersections, deallocate
  use libsupermesh_intersection_finder, only : &
    & advancing_front_intersection_finder
  use libsupermesh_read_triangle, only : read_ele, read_node
  use libsupermesh_unittest_tools, only : report_test

  implicit none
  
  integer :: i
  logical :: fail
  integer, dimension(:, :), allocatable :: enlist_a, enlist_b
  real, dimension(:, :), allocatable :: positions_a, positions_b
  type(intersections), dimension(1) :: map_ab

  integer, parameter :: dim = 3

  call read_node("data/tet.node", dim, positions_a)
  call read_ele("data/tet.ele", dim, enlist_a)
  call read_node("data/tet.node", dim, positions_b)
  call read_ele("data/tet.ele", dim, enlist_b)
  call advancing_front_intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)

  fail = (map_ab(1)%n /= 1)
  call report_test("[intersection finder: length]", fail, .false., "There shall be only one element")

  i = map_ab(1)%v(1)
  fail = (i /= 1)
  call report_test("[intersection finder: correct]", fail, .false., "The answer should be one")
    
  call deallocate(map_ab)
  deallocate(positions_a, enlist_a, positions_b, enlist_b)

end subroutine test_intersection_finder_3d
