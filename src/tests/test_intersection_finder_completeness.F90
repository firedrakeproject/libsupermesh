subroutine test_intersection_finder_completeness

  use libsupermesh_intersection_finder_module
  use libsupermesh_fields_dummy
  use libsupermesh_read_triangle_2
  use libsupermesh_unittest_tools
  use libsupermesh_linked_lists
  use libsupermesh_construction
  use libsupermesh_tri_intersection_module
  use libsupermesh_tet_intersection_module, only : tet_buf_size

  implicit none

  type(vector_field) :: positionsA, positionsB
  type(ilist), dimension(:), allocatable :: map_BA
  integer :: ele_A, ele_B, ele_C, n_trisC, n_tetsC, i
  real :: vol_B, vols_C
  logical :: fail
  type(inode), pointer :: llnode
  type(vector_field) :: intersection
  real, dimension(2, 3, tri_buf_size) ::  trisC_real
  real, dimension(3, 4, tet_buf_size) ::  tetsC_real
  type(mesh_type) :: new_mesh, intersection_mesh

  integer, parameter :: dim = 2, loc = 3

  positionsA = read_triangle_files("data/intersection_finder.1", dim)
  positionsB = read_triangle_files("data/intersection_finder.2", dim)

  allocate(map_BA(ele_count(positionsB)))

  map_BA = advancing_front_intersection_finder( &
      & positionsB%val, reshape(positionsB%mesh%ndglno, (/loc, ele_count(positionsB)/)), &
      & positionsA%val, reshape(positionsA%mesh%ndglno, (/loc, ele_count(positionsA)/)) )
  call libsupermesh_cintersector_set_dimension(dim)

  do ele_B=1,ele_count(positionsB)
    vol_B = triangle_area(ele_val(positionsB, ele_B))

    llnode => map_BA(ele_B)%firstnode
    vols_C = 0.0
    do while(associated(llnode))
      ele_A = llnode%value

      ! Triangles (2D)
        call intersect_elements(ele_val(positionsA, ele_A), ele_val(positionsB, ele_B), &
                    n_trisC, trisC_real)
        call allocate(new_mesh, dim, n_trisC * loc, n_trisC, loc)

        if ( n_trisC > 0 ) then
          new_mesh%ndglno = (/ (i, i=1,loc * n_trisC) /)
          new_mesh%continuity = -1
        end if

        call allocate(intersection, dim, new_mesh)
        if ( n_trisC > 0 ) then
          do i = 1, n_trisC
            call set(intersection, ele_nodes(intersection, i), trisC_real(:,:,i))
          end do
        end if
        call deallocate(new_mesh)

      do ele_C=1,ele_count(intersection)
        vols_C = vols_C + triangle_area(ele_val(intersection, ele_C))
      end do
      llnode => llnode%next
    end do

    fail = (vol_B .fne. vols_C)
    call report_test("[intersection finder: completeness]", fail, .false., "Need to have the same volume!")
    if (fail) then
      write(0,*) "ele_B: ", ele_B
      write(0,*) "vol_B: ", vol_B
      write(0,*) "vols_C: ", vols_C
    end if
    call deallocate(intersection)
  end do
    
end subroutine test_intersection_finder_completeness
