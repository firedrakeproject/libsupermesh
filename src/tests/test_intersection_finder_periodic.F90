subroutine test_intersection_finder_periodic

  use libsupermesh_intersection_finder
  use libsupermesh_fields_dummy
  use libsupermesh_read_triangle_2
  use libsupermesh_unittest_tools
  use libsupermesh_linked_lists
  implicit none

  type(vector_field) :: posA, posB
  type(ilist), dimension(1) :: map_AB
  logical :: fail

  integer, parameter :: dim = 2
  integer :: i = 0

  ! A has one element
  ! B has two disconnected elements

  posA = read_triangle_files("data/intersection_finder_periodic_A", dim)
  posB = read_triangle_files("data/intersection_finder_periodic_B", dim)

  call intersection_finder( &
      & posA%val, reshape(posA%mesh%ndglno, (/ele_loc(posA, 1), ele_count(posA)/)), &
      & posB%val, reshape(posB%mesh%ndglno, (/ele_loc(posB, 1), ele_count(posB)/)), map_AB)

  fail = (map_AB(1)%length /= 2)
  call report_test("[intersection finder periodic]", fail, .false., "")
  
  call deallocate(posA)
  call deallocate(posB)
  
  do i=1,size(map_AB)
    call flush_list(map_AB(i))
  end do

end subroutine test_intersection_finder_periodic
