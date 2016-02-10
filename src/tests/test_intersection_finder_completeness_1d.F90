#include "libsupermesh_debug.h"

subroutine test_intersection_finder_completeness_1d() bind(c)

  use libsupermesh_intersection_finder, only : intersections, deallocate, &
    & intersection_finder, advancing_front_intersection_finder, &
    & sort_intersection_finder, brute_force_intersection_finder
  use libsupermesh_interval_intersection, only : interval_buf_size, &
    & intersect_intervals
  use libsupermesh_read_triangle, only : read_ele, read_node
  use libsupermesh_unittest_tools, only : report_test, operator(.fne.)

  implicit none

  integer :: ele_a, ele_b, ele_c, i, n_intervalsC
  real :: size_b, size_c
  integer, dimension(:, :), allocatable :: enlist_a, enlist_b
  logical :: fail
  real, dimension(1, 2, interval_buf_size) :: intervalsC
  real, dimension(:, :), allocatable :: positions_a, positions_b
  type(intersections), dimension(:), allocatable :: map_ba

  integer, parameter :: dim = 1

  call read_node("data/line.1.node", dim, positions_a)
  call read_ele("data/line.1.ele", dim, enlist_a)
  call read_node("data/line.2.node", dim, positions_b)
  call read_ele("data/line.2.ele", dim, enlist_b)
  allocate(map_ba(size(enlist_b, 2)))
  
  call intersection_finder(positions_b, enlist_b, positions_a, enlist_a, map_ba)
  call check_map_ba("intersection_finder")
  call deallocate(map_ba)
  
  call advancing_front_intersection_finder(positions_b, enlist_b, positions_a, enlist_a, map_ba)
  call check_map_ba("advancing_front_intersection_finder")
  call deallocate(map_ba)
  
  call sort_intersection_finder(positions_b, enlist_b, positions_a, enlist_a, map_ba)
  call check_map_ba("sort_intersection_finder")
  call deallocate(map_ba)
  
  call brute_force_intersection_finder(positions_b, enlist_b, positions_a, enlist_a, map_ba)
  call check_map_ba("brute_force_intersection_finder")
  call deallocate(map_ba)

  deallocate(map_ba, positions_a, enlist_a, positions_b, enlist_b)

contains

  pure function interval_size(interval) result(size)
    real, dimension(1, 2), intent(in) :: interval
    
    real :: size
    
    size = abs(interval(1, 2) - interval(1, 1))
    
  end function interval_size

  subroutine check_map_ba(name)
    character(len = *), intent(in) :: name

    fail = .false.
    do ele_b = 1, size(enlist_b, 2)
      size_b = interval_size(positions_b(:, enlist_b(:, ele_b)))

      size_c = 0.0
      do i = 1, map_ba(ele_b)%n
        ele_a = map_ba(ele_b)%v(i)
        call intersect_intervals(positions_a(:, enlist_a(:, ele_a)), positions_b(:, enlist_b(:, ele_b)), intervalsC, n_intervalsC)
        do ele_c = 1, n_intervalsC
          size_c = size_c + interval_size(intervalsC(:, :, ele_c))
        end do
      end do

      fail = (size_b .fne. size_c)
      if(fail) exit
    end do
    call report_test("[" // trim(name) // ": completeness]", fail, .false., "Need to have the same size")
 
 end subroutine check_map_ba

end subroutine test_intersection_finder_completeness_1d
