#include "libsupermesh_debug.h"

subroutine test_tet_intersector() bind(c)

  use libsupermesh_intersection_finder, only : intersections, deallocate, &
    & intersection_finder
  use libsupermesh_supermesh, only : intersect_elements, &
    & intersect_simplices, tetrahedron_volume
  use libsupermesh_read_triangle, only : read_ele, read_node
  use libsupermesh_tet_intersection, only : tet_type, tet_buf_size, &
    & intersect_tets, get_planes
  use libsupermesh_unittest_tools, only : report_test, operator(.fne.)
  
  implicit none

  integer :: ele_a, ele_b, ele_c, i, n_tets_c
  integer, dimension(:, :), allocatable :: enlist_a, enlist_b
  real :: volume_c
  real, dimension(3, 4) :: tet_a_real
  real, dimension(3, 4, tet_buf_size) :: tets_c_real
  real, dimension(:, :), allocatable :: positions_a, positions_b 
  type(intersections), dimension(:), allocatable :: map_ab
  type(tet_type) :: tet_a, tet_b
  type(tet_type), dimension(tet_buf_size) :: tets_c, work

  integer, parameter :: dim = 3
  
  call read_node("data/pyramid_0_05.node", dim, positions_a)
  call read_ele("data/pyramid_0_05.ele", dim, enlist_a)
  call read_node("data/cube_0_05.node", dim, positions_b)
  call read_ele("data/cube_0_05.ele", dim, enlist_b)
  
  allocate(map_ab(size(enlist_a, 2)))
  call intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)

  volume_c = 0.0D0
  do ele_a = 1, size(enlist_a, 2)
    tet_a%v = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      tet_b%v = positions_b(:, enlist_b(:, ele_b))
      call intersect_tets(tet_a, tet_b, tets_c, n_tets_c)
      do ele_c = 1, n_tets_c
        volume_c = volume_c + tetrahedron_volume(tets_c(ele_c)%v)
      end do
    end do    
  end do
  call report_test("[intersect_tets]", volume_c .fne. 1.0D0 / 3.0D0, .false., "Incorrect intersection volume")

  volume_c = 0.0D0
  do ele_a = 1, size(enlist_a, 2)
    tet_a_real = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_tets(tet_a_real, positions_b(:, enlist_b(:, ele_b)), tets_c_real, n_tets_c)
      do ele_c = 1, n_tets_c
        volume_c = volume_c + tetrahedron_volume(tets_c_real(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_tets]", volume_c .fne. 1.0D0 / 3.0D0, .false., "Incorrect intersection volume")

  volume_c = 0.0D0
  do ele_a = 1, size(enlist_a, 2)
    tet_a%v = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      tet_b%v = positions_b(:, enlist_b(:, ele_b))
      call intersect_tets(tet_a, get_planes(tet_b), tets_c, n_tets_c, work = work)
      do ele_c = 1, n_tets_c
        volume_c = volume_c + tetrahedron_volume(tets_c(ele_c)%v)
      end do
    end do    
  end do
  call report_test("[intersect_tets]", volume_c .fne. 1.0D0 / 3.0D0, .false., "Incorrect intersection volume")
  
  volume_c = 0.0D0
  do ele_a = 1, size(enlist_a, 2)
    tet_a_real = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_simplices(tet_a_real, positions_b(:, enlist_b(:, ele_b)), tets_c_real, n_tets_c)
      do ele_c = 1, n_tets_c
        volume_c = volume_c + tetrahedron_volume(tets_c_real(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_simplices]", volume_c .fne. 1.0D0 / 3.0D0, .false., "Incorrect intersection volume")
  
  volume_c = 0.0D0
  do ele_a = 1, size(enlist_a, 2)
    tet_a_real = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_elements(tet_a_real, positions_b(:, enlist_b(:, ele_b)), tets_c_real, n_tets_c)
      do ele_c = 1, n_tets_c
        volume_c = volume_c + tetrahedron_volume(tets_c_real(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_elements]", volume_c .fne. 1.0D0 / 3.0D0, .false., "Incorrect intersection volume")
  
  call deallocate(map_ab)
  deallocate(map_ab, positions_a, enlist_a, positions_b, enlist_b)

end subroutine test_tet_intersector
