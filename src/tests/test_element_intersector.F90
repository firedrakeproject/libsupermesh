#include "libsupermesh_debug.h"

subroutine test_element_intersector() bind(c)

  use libsupermesh_intersection_finder, only : intersections, deallocate, &
    & intersection_finder
  use libsupermesh_supermesh, only : max_n_simplices_c, intersect_elements, &
    & simplex_volume
  use libsupermesh_read_triangle, only : read_ele, read_node
  use libsupermesh_unittest_tools, only : report_test, operator(.fne.)
  
  implicit none

  integer :: ele_a, ele_b, ele_c, i, n_elementsC
  integer, dimension(:, :), allocatable :: enlist_a, enlist_b
  real :: volume_c
  real, dimension(:, :), allocatable :: elementA, positions_a, positions_b
  real, dimension(:, :, :), allocatable :: elementsC
  type(intersections), dimension(:), allocatable :: map_ab
  
  ! Interval-interval
  
  call read_node("data/line.1.node", 1, positions_a)
  call read_ele("data/line.1.ele", 1, enlist_a)
  call read_node("data/line.2.node", 1, positions_b)
  call read_ele("data/line.2.ele", 1, enlist_b)
  allocate(elementA(size(positions_a, 1), size(enlist_a, 1)), &
    & elementsC(size(positions_a, 1), size(positions_a, 1) + 1, &
              & max_n_simplices_c(size(positions_a, 1), size(enlist_a, 1), size(enlist_b, 1))))
  
  allocate(map_ab(size(enlist_a, 2)))
  call intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
  
  volume_c = 0.0D0
  do ele_a = 1, size(enlist_a, 2)
    elementA = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_elements(elementA, positions_b(:, enlist_b(:, ele_b)), elementsC, n_elementsC)
      do ele_c = 1, n_elementsC
        volume_c = volume_c + simplex_volume(elementsC(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_elements]", volume_c .fne. 1.0D0, .false., "Incorrect intersection volume")
  
  call deallocate(map_ab)
  deallocate(map_ab, positions_a, enlist_a, positions_b, enlist_b, elementA, elementsC)
  
  ! Tri-tri
  
  call read_node("data/square.1.node", 2, positions_a)
  call read_ele("data/square.1.ele", 2, enlist_a)
  call read_node("data/square.2.node", 2, positions_b)
  call read_ele("data/square.2.ele", 2, enlist_b)
  allocate(elementA(size(positions_a, 1), size(enlist_a, 1)), &
    & elementsC(size(positions_a, 1), size(positions_a, 1) + 1, &
              & max_n_simplices_c(size(positions_a, 1), size(enlist_a, 1), size(enlist_b, 1))))
  
  allocate(map_ab(size(enlist_a, 2)))
  call intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
  
  volume_c = 0.0D0
  do ele_a = 1, size(enlist_a, 2)
    elementA = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_elements(elementA, positions_b(:, enlist_b(:, ele_b)), elementsC, n_elementsC)
      do ele_c = 1, n_elementsC
        volume_c = volume_c + simplex_volume(elementsC(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_elements]", volume_c .fne. 1.0D0, .false., "Incorrect intersection volume")
  
  call deallocate(map_ab)
  deallocate(map_ab, positions_a, enlist_a, positions_b, enlist_b, elementA, elementsC)
  
  ! Tri-quad
  
  call read_node("data/square.1.node", 2, positions_a)
  call read_ele("data/square.1.ele", 2, enlist_a)
  call read_node("data/square.4.node", 2, positions_b)
  call read_ele("data/square.4.ele", 2, enlist_b)
  allocate(elementA(size(positions_a, 1), size(enlist_a, 1)), &
    & elementsC(size(positions_a, 1), size(positions_a, 1) + 1, &
              & max_n_simplices_c(size(positions_a, 1), size(enlist_a, 1), size(enlist_b, 1))))
  
  allocate(map_ab(size(enlist_a, 2)))
  call intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
  
  volume_c = 0.0D0
  do ele_a = 1, size(enlist_a, 2)
    elementA = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_elements(elementA, positions_b(:, enlist_b(:, ele_b)), elementsC, n_elementsC)
      do ele_c = 1, n_elementsC
        volume_c = volume_c + simplex_volume(elementsC(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_elements]", volume_c .fne. 1.0D0, .false., "Incorrect intersection volume")
  
  call deallocate(map_ab)
  deallocate(map_ab, positions_a, enlist_a, positions_b, enlist_b, elementA, elementsC)
  
  ! Quad-tri
  
  call read_node("data/square.3.node", 2, positions_a)
  call read_ele("data/square.3.ele", 2, enlist_a)
  call read_node("data/square.2.node", 2, positions_b)
  call read_ele("data/square.2.ele", 2, enlist_b)
  allocate(elementA(size(positions_a, 1), size(enlist_a, 1)), &
    & elementsC(size(positions_a, 1), size(positions_a, 1) + 1, &
              & max_n_simplices_c(size(positions_a, 1), size(enlist_a, 1), size(enlist_b, 1))))
  
  allocate(map_ab(size(enlist_a, 2)))
  call intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
  
  volume_c = 0.0D0
  do ele_a = 1, size(enlist_a, 2)
    elementA = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_elements(elementA, positions_b(:, enlist_b(:, ele_b)), elementsC, n_elementsC)
      do ele_c = 1, n_elementsC
        volume_c = volume_c + simplex_volume(elementsC(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_elements]", volume_c .fne. 1.0D0, .false., "Incorrect intersection volume")
  
  call deallocate(map_ab)
  deallocate(map_ab, positions_a, enlist_a, positions_b, enlist_b, elementA, elementsC)
  
  ! Quad-quad
  
  call read_node("data/square.3.node", 2, positions_a)
  call read_ele("data/square.3.ele", 2, enlist_a)
  call read_node("data/square.4.node", 2, positions_b)
  call read_ele("data/square.4.ele", 2, enlist_b)
  allocate(elementA(size(positions_a, 1), size(enlist_a, 1)), &
    & elementsC(size(positions_a, 1), size(positions_a, 1) + 1, &
              & max_n_simplices_c(size(positions_a, 1), size(enlist_a, 1), size(enlist_b, 1))))
  
  allocate(map_ab(size(enlist_a, 2)))
  call intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
  
  volume_c = 0.0D0
  do ele_a = 1, size(enlist_a, 2)
    elementA = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_elements(elementA, positions_b(:, enlist_b(:, ele_b)), elementsC, n_elementsC)
      do ele_c = 1, n_elementsC
        volume_c = volume_c + simplex_volume(elementsC(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_elements]", volume_c .fne. 1.0D0, .false., "Incorrect intersection volume")
  
  call deallocate(map_ab)
  deallocate(map_ab, positions_a, enlist_a, positions_b, enlist_b, elementA, elementsC)
  
  ! Tet-tet
  
  call read_node("data/cube.1.node", 3, positions_a)
  call read_ele("data/cube.1.ele", 3, enlist_a)
  call read_node("data/cube.2.node", 3, positions_b)
  call read_ele("data/cube.2.ele", 3, enlist_b)
  allocate(elementA(size(positions_a, 1), size(enlist_a, 1)), &
    & elementsC(size(positions_a, 1), size(positions_a, 1) + 1, &
              & max_n_simplices_c(size(positions_a, 1), size(enlist_a, 1), size(enlist_b, 1))))
  
  allocate(map_ab(size(enlist_a, 2)))
  call intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
  
  volume_c = 0.0D0
  do ele_a = 1, size(enlist_a, 2)
    elementA = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_elements(elementA, positions_b(:, enlist_b(:, ele_b)), elementsC, n_elementsC)
      do ele_c = 1, n_elementsC
        volume_c = volume_c + simplex_volume(elementsC(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_elements]", volume_c .fne. 1.0D0, .false., "Incorrect intersection volume")
  
  call deallocate(map_ab)
  deallocate(map_ab, positions_a, enlist_a, positions_b, enlist_b, elementA, elementsC)
  
  ! Tet-hex
  
  call read_node("data/cube.1.node", 3, positions_a)
  call read_ele("data/cube.1.ele", 3, enlist_a)
  call read_node("data/cube.3.node", 3, positions_b)
  call read_ele("data/cube.3.ele", 3, enlist_b)
  allocate(elementA(size(positions_a, 1), size(enlist_a, 1)), &
    & elementsC(size(positions_a, 1), size(positions_a, 1) + 1, &
              & max_n_simplices_c(size(positions_a, 1), size(enlist_a, 1), size(enlist_b, 1))))
  
  allocate(map_ab(size(enlist_a, 2)))
  call intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
  
  volume_c = 0.0D0
  do ele_a = 1, size(enlist_a, 2)
    elementA = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_elements(elementA, positions_b(:, enlist_b(:, ele_b)), elementsC, n_elementsC)
      do ele_c = 1, n_elementsC
        volume_c = volume_c + simplex_volume(elementsC(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_elements]", volume_c .fne. 0.5D0, .false., "Incorrect intersection volume")
  
  call deallocate(map_ab)
  deallocate(map_ab, positions_a, enlist_a, positions_b, enlist_b, elementA, elementsC)
  
  ! Hex-tet
  
  call read_node("data/cube.3.node", 3, positions_a)
  call read_ele("data/cube.3.ele", 3, enlist_a)
  call read_node("data/cube.2.node", 3, positions_b)
  call read_ele("data/cube.2.ele", 3, enlist_b)
  allocate(elementA(size(positions_a, 1), size(enlist_a, 1)), &
    & elementsC(size(positions_a, 1), size(positions_a, 1) + 1, &
              & max_n_simplices_c(size(positions_a, 1), size(enlist_a, 1), size(enlist_b, 1))))
  
  allocate(map_ab(size(enlist_a, 2)))
  call intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
  
  volume_c = 0.0D0
  do ele_a = 1, size(enlist_a, 2)
    elementA = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_elements(elementA, positions_b(:, enlist_b(:, ele_b)), elementsC, n_elementsC)
      do ele_c = 1, n_elementsC
        volume_c = volume_c + simplex_volume(elementsC(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_elements]", volume_c .fne. 0.5D0, .false., "Incorrect intersection volume")
  
  call deallocate(map_ab)
  deallocate(map_ab, positions_a, enlist_a, positions_b, enlist_b, elementA, elementsC)
  
  ! Hex-hex
  
  call read_node("data/cube.3.node", 3, positions_a)
  call read_ele("data/cube.3.ele", 3, enlist_a)
  call read_node("data/cube.4.node", 3, positions_b)
  call read_ele("data/cube.4.ele", 3, enlist_b)
  allocate(elementA(size(positions_a, 1), size(enlist_a, 1)), &
    & elementsC(size(positions_a, 1), size(positions_a, 1) + 1, &
              & max_n_simplices_c(size(positions_a, 1), size(enlist_a, 1), size(enlist_b, 1))))
  
  allocate(map_ab(size(enlist_a, 2)))
  call intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
  
  volume_c = 0.0D0
  do ele_a = 1, size(enlist_a, 2)
    elementA = positions_a(:, enlist_a(:, ele_a))
    do i = 1, map_ab(ele_a)%n
      ele_b = map_ab(ele_a)%v(i)
      call intersect_elements(elementA, positions_b(:, enlist_b(:, ele_b)), elementsC, n_elementsC)
      do ele_c = 1, n_elementsC
        volume_c = volume_c + simplex_volume(elementsC(:, :, ele_c))
      end do
    end do    
  end do
  call report_test("[intersect_elements]", volume_c .fne. 0.5D0, .false., "Incorrect intersection volume")
  
  call deallocate(map_ab)
  deallocate(map_ab, positions_a, enlist_a, positions_b, enlist_b, elementA, elementsC)

end subroutine test_element_intersector
