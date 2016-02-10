#include "libsupermesh_debug.h"

subroutine test_tet_intersector()

  use libsupermesh_intersection_finder, only : intersections, deallocate, &
    & intersection_finder
  use libsupermesh_supermesh, only : intersect_elements, &
    & intersect_simplices, tetrahedron_volume
  use libsupermesh_read_triangle, only : read_ele, read_node
  use libsupermesh_tet_intersection, only : tet_type, tet_buf_size, &
    & intersect_tets, get_planes
  use libsupermesh_unittest_tools, only : report_test, operator(.fne.)
  
  implicit none

  integer :: ele_a, ele_b, ele_c, i, n_tetsC
  integer, dimension(:, :), allocatable :: enlist_a, enlist_b
  real :: volume_c
  real, dimension(3, 4) :: tetA_real
  real, dimension(3, 4, tet_buf_size) :: tetsC_real
  real, dimension(:, :), allocatable :: positions_a, positions_b 
  type(intersections), dimension(:), allocatable :: map_ab
  type(tet_type) :: tetA, tetB
  type(tet_type), dimension(tet_buf_size) :: tetsC, work

  integer, parameter :: dim = 3
  
  call read_node("data/pyramid_0_9_4.node", dim, positions_a)
  call read_ele("data/pyramid_0_9_4.ele", dim, enlist_a)
  call read_node("data/cube_0_9_4.node", dim, positions_b)
  call read_ele("data/cube_0_9_4.ele", dim, enlist_b)
  
  allocate(map_ab(size(enlist_a, 2)))
  call intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)

  volume_c = 0.0
  do ele_a = 1, size(enlist_a, 2)
    tetA%v = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      tetB%v = positions_b(:, enlist_b(:, ele_b))
      call intersect_tets(tetA, tetB, tetsC, n_tetsC)
      do ele_c = 1, n_tetsC
        volume_c = volume_c + tetrahedron_volume(tetsC(ele_c)%v)
      end do
    end do    
  end do
  call report_test("[intersect_tets]", volume_c .fne. 1000.0 / 3.0, .false., "Incorrect intersection volume")

  volume_c = 0.0
  do ele_a = 1, size(enlist_a, 2)
    tetA_real = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_tets(tetA_real, positions_b(:, enlist_b(:, ele_b)), tetsC_real, n_tetsC)
      do ele_c = 1, n_tetsC
        volume_c = volume_c + tetrahedron_volume(tetsC_real(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_tets]", volume_c .fne. 1000.0 / 3.0, .false., "Incorrect intersection volume")

  volume_c = 0.0
  do ele_a = 1, size(enlist_a, 2)
    tetA%v = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      tetB%v = positions_b(:, enlist_b(:, ele_b))
      call intersect_tets(tetA, get_planes(tetB), tetsC, n_tetsC, work = work)
      do ele_c = 1, n_tetsC
        volume_c = volume_c + tetrahedron_volume(tetsC(ele_c)%v)
      end do
    end do    
  end do
  call report_test("[intersect_tets]", volume_c .fne. 1000.0 / 3.0, .false., "Incorrect intersection volume")
  
  volume_c = 0.0
  do ele_a = 1, size(enlist_a, 2)
    tetA_real = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_simplices(tetA_real, positions_b(:, enlist_b(:, ele_b)), tetsC_real, n_tetsC)
      do ele_c = 1, n_tetsC
        volume_c = volume_c + tetrahedron_volume(tetsC_real(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_simplices]", volume_c .fne. 1000.0 / 3.0, .false., "Incorrect intersection volume")
  
  volume_c = 0.0
  do ele_a = 1, size(enlist_a, 2)
    tetA_real = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_elements(tetA_real, positions_b(:, enlist_b(:, ele_b)), tetsC_real, n_tetsC)
      do ele_c = 1, n_tetsC
        volume_c = volume_c + tetrahedron_volume(tetsC_real(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_elements]", volume_c .fne. 1000.0 / 3.0, .false., "Incorrect intersection volume")
  
  call deallocate(map_ab)
  deallocate(map_ab, positions_a, enlist_a, positions_b, enlist_b)

end subroutine test_tet_intersector
