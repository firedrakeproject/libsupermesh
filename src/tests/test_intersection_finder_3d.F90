subroutine test_intersection_finder_3d

  use libsupermesh_intersection_finder_module
  use libsupermesh_fields
  use libsupermesh_read_triangle
  use libsupermesh_unittest_tools
  use libsupermesh_linked_lists
  
  type(vector_field) :: positionsA, positionsB
  type(ilist), dimension(1) :: map_AB

  integer :: i
  logical :: fail

  positionsA = read_triangle_files("data/tet", quad_degree=4)
  positionsB = read_triangle_files("data/tet", quad_degree=4)
  map_AB = libsupermesh_advancing_front_intersection_finder( &
      & positionsA%val, reshape(positionsA%mesh%ndglno, (/ele_loc(positionsA, 1), ele_count(positionsA)/)), &
      & positionsB%val, reshape(positionsB%mesh%ndglno, (/ele_loc(positionsB, 1), ele_count(positionsB)/)) )

  fail = (map_AB(1)%length /= 1)
  call report_test("[intersection finder: length]", fail, .false., "There shall be only one")

  i = fetch(map_AB(1), 1)
  fail = (i /= 1)
  call report_test("[intersection finder: correct]", fail, .false., "The answer should be one")

end subroutine test_intersection_finder_3d
