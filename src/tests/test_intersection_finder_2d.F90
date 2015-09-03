subroutine test_intersection_finder_2d

  use libsupermesh_intersection_finder
  use libsupermesh_fields
  use libsupermesh_read_triangle
  use libsupermesh_unittest_tools
  use libsupermesh_linked_lists
  
  implicit none

  type(vector_field) :: positionsA, positionsB
  type(ilist), dimension(1) :: map_AB
  type(ilist), dimension(3) :: bigger_map_AB
  !type(inode), pointer :: node

  integer :: i
  logical :: fail

  integer, parameter :: dim = 2, loc = 3
  
  positionsA = read_triangle_files("data/triangle.1", dim)
  positionsB = read_triangle_files("data/triangle.1", dim)
  call intersection_finder( &
      & positionsA%val, reshape(positionsA%mesh%ndglno, (/loc, ele_count(positionsA)/)), &
      & positionsB%val, reshape(positionsB%mesh%ndglno, (/loc, ele_count(positionsB)/)), map_AB)

  fail = (map_AB(1)%length /= 1)
  call report_test("[intersection finder: length]", fail, .false., "There shall be only one")

  i = fetch(map_AB(1), 1)
  fail = (i /= 1)
  call report_test("[intersection finder: correct]", fail, .false., "The answer should be one")
  
  do i=1,size(map_AB)
    call flush_list(map_AB(i))
  end do

  call deallocate(positionsB)
  positionsB = read_triangle_files("data/triangle.2", dim)
  call intersection_finder( &
      & positionsA%val, reshape(positionsA%mesh%ndglno, (/loc, ele_count(positionsA)/)), &
      & positionsB%val, reshape(positionsB%mesh%ndglno, (/loc, ele_count(positionsB)/)), map_AB)

  fail = (map_AB(1)%length /= 3)
  call report_test("[intersection finder: length]", fail, .false., "There shall be three elements")

  call deallocate(positionsA)
  positionsA = read_triangle_files("data/triangle.2", dim)
  call intersection_finder( &
      & positionsA%val, reshape(positionsA%mesh%ndglno, (/loc, ele_count(positionsA)/)), &
      & positionsB%val, reshape(positionsB%mesh%ndglno, (/loc, ele_count(positionsB)/)), bigger_map_AB)
  do i=1,ele_count(positionsA)
    fail = (bigger_map_AB(i)%length < 1)
    call report_test("[intersection finder: length]", fail, .false., "There shall be only one")

    fail = (.not. has_value(bigger_map_AB(i), i))
    call report_test("[intersection finder: correct]", fail, .false., "The answer should be correct")
  end do

  call deallocate(positionsA)
  call deallocate(positionsB)

  do i=1,size(map_AB)
    call flush_list(map_AB(i))
  end do
  
  do i=1,size(bigger_map_AB)
    call flush_list(bigger_map_AB(i))
  end do

end subroutine test_intersection_finder_2d
