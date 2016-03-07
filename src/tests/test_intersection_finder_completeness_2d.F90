#include "libsupermesh_debug.h"

subroutine test_intersection_finder_completeness_2d() bind(c)

  use libsupermesh_intersection_finder, only : intersections, deallocate, &
    & intersection_finder, advancing_front_intersection_finder, &
    & quadtree_intersection_finder, tree_intersection_finder, &
    & rtree_intersection_finder, brute_force_intersection_finder
  use libsupermesh_supermesh, only : triangle_area
  use libsupermesh_read_triangle, only : read_ele, read_node
  use libsupermesh_tri_intersection, only : intersect_tris, tri_buf_size
  use libsupermesh_unittest_tools, only : report_test, operator(.fne.)

  implicit none

  integer :: ele_a, ele_b, ele_c, i, n_tris_c
  real :: area_b, area_c
  integer, dimension(:, :), allocatable :: enlist_a, enlist_b
  logical :: fail
  real, dimension(2, 3, tri_buf_size) ::  tris_c
  real, dimension(:, :), allocatable :: positions_a, positions_b
  type(intersections), dimension(:), allocatable :: map_ba

  integer, parameter :: dim = 2

  call read_node("data/square.1.node", dim, positions_a)
  call read_ele("data/square.1.ele", dim, enlist_a)
  call read_node("data/square.2.node", dim, positions_b)
  call read_ele("data/square.2.ele", dim, enlist_b)
  allocate(map_ba(size(enlist_b, 2)))
  
  call intersection_finder(positions_b, enlist_b, positions_a, enlist_a, map_ba)
  call check_map_ba("intersection_finder")
  call deallocate(map_ba)
  
  call advancing_front_intersection_finder(positions_b, enlist_b, positions_a, enlist_a, map_ba)
  call check_map_ba("advancing_front_intersection_finder")
  call deallocate(map_ba)
  
  call quadtree_intersection_finder(positions_b, enlist_b, positions_a, enlist_a, map_ba)
  call check_map_ba("quadtree_intersection_finder")
  call deallocate(map_ba)
  
  call rtree_intersection_finder(positions_b, enlist_b, positions_a, enlist_a, map_ba)
  call check_map_ba("rtree_intersection_finder")
  call deallocate(map_ba)
  
  call tree_intersection_finder(positions_b, enlist_b, positions_a, enlist_a, map_ba)
  call check_map_ba("tree_intersection_finder")
  call deallocate(map_ba)
  
  call brute_force_intersection_finder(positions_b, enlist_b, positions_a, enlist_a, map_ba)
  call check_map_ba("brute_force_intersection_finder")
  call deallocate(map_ba)

  deallocate(map_ba, positions_a, enlist_a, positions_b, enlist_b)

contains

  subroutine check_map_ba(name)
    character(len = *), intent(in) :: name

    fail = .false.
    do ele_b = 1, size(enlist_b, 2)
      area_b = triangle_area(positions_b(:, enlist_b(:, ele_b)))

      area_c = 0.0D0
      do i = 1, map_ba(ele_b)%n
        ele_a = map_ba(ele_b)%v(i)
        call intersect_tris(positions_a(:, enlist_a(:, ele_a)), positions_b(:, enlist_b(:, ele_b)), tris_c, n_tris_c)
        do ele_c = 1, n_tris_c
          area_c = area_c + triangle_area(tris_c(:, :, ele_c))
        end do
      end do

      fail = (area_b .fne. area_c)
      if(fail) exit
    end do
    call report_test("[" // trim(name) // ": completeness]", fail, .false., "Need to have the same area")
 
 end subroutine check_map_ba

end subroutine test_intersection_finder_completeness_2d
