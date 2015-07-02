subroutine test_intersection_finder_3d

  use libsupermesh_intersection_finder_module
  use libsupermesh_fields_dummy
  use libsupermesh_read_triangle_2
  use libsupermesh_unittest_tools
  use libsupermesh_linked_lists

  implicit none
  
  type(vector_field) :: positionsA, positionsB
  type(ilist), dimension(1) :: map_AB

  integer :: i
  logical :: fail

  integer, parameter :: dim = 3, loc = 4

  positionsA = read_triangle_files("data/tet", dim)
  positionsB = read_triangle_files("data/tet", dim)
  map_AB = advancing_front_intersection_finder( &
      & positionsA%val, reshape(positionsA%mesh%ndglno, (/loc, ele_count(positionsA)/)), &
      & positionsB%val, reshape(positionsB%mesh%ndglno, (/loc, ele_count(positionsB)/)) )

  fail = (map_AB(1)%length /= 1)
  call report_test("[intersection finder: length]", fail, .false., "There shall be only one")

  i = fetch(map_AB(1), 1)
  fail = (i /= 1)
  call report_test("[intersection finder: correct]", fail, .false., "The answer should be one")

end subroutine test_intersection_finder_3d
