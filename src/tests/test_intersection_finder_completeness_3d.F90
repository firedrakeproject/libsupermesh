#include "confdefs.h"

subroutine test_intersection_finder_completeness_3d
  
  use libsupermesh_intersection_finder_module
  use libsupermesh_fields_dummy
  use libsupermesh_read_triangle_2
  use libsupermesh_unittest_tools
  use libsupermesh_linked_lists
  use libsupermesh_elements
  use libsupermesh_construction
  use libsupermesh_tri_intersection_module, only : tri_buf_size
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

  integer, parameter :: dim = 3, loc = 4

  positionsA = read_triangle_files("data/cube.1", dim)
  positionsB = read_triangle_files("data/cube.2", dim)

  allocate(map_BA(ele_count(positionsB)))

  map_BA = advancing_front_intersection_finder( &
      & positionsB%val, reshape(positionsB%mesh%ndglno, (/loc, ele_count(positionsB)/)), &
      & positionsA%val, reshape(positionsA%mesh%ndglno, (/loc, ele_count(positionsA)/)) )
  call intersector_set_dimension(dim)

  do ele_B=1,ele_count(positionsB)
    vol_B = tetrahedron_volume(ele_val(positionsB, ele_B))

    llnode => map_BA(ele_B)%firstnode
    vols_C = 0.0
    do while(associated(llnode))
      ele_A = llnode%value
      ! Tets (3D)
        call intersect_elements(ele_val(positionsA, ele_A), ele_val(positionsB, ele_B), &
                dim, n_tetsC, tetsC_real=tetsC_real )
        call allocate(intersection_mesh, n_tetsC * loc, n_tetsC, loc)
        intersection_mesh%continuity = -1

        if ( n_tetsC > 0 ) then
          intersection_mesh%ndglno = (/ (i, i=1,loc * n_tetsC) /)
        end if
 
        call allocate(intersection, dim, intersection_mesh)
        if ( n_tetsC > 0 ) then
          do i = 1, n_tetsC
            call set(intersection, ele_nodes(intersection, i), tetsC_real(:,:,i))
          end do
        end if
        call deallocate(intersection_mesh)

      do ele_C=1,ele_count(intersection)
        vols_C = vols_C + tetrahedron_volume(ele_val(intersection, ele_C))
      end do
      llnode => llnode%next
    end do

    fail = (vol_B .fne. vols_C)
    call report_test("[intersection finder 3D: completeness]", fail, .false., "Need to have the same volume!")
    if (fail) then
      write(0,*) "ele_B: ", ele_B
      write(0,*) "vol_B: ", vol_B
      write(0,*) "vols_C: ", vols_C
    end if
  end do
    
end subroutine test_intersection_finder_completeness_3d
