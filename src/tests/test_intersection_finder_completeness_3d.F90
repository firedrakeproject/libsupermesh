#include "confdefs.h"

subroutine test_intersection_finder_completeness_3d
  
  use libsupermesh_intersection_finder_module
  use libsupermesh_fields
  use libsupermesh_read_triangle
  use libsupermesh_unittest_tools
  use libsupermesh_linked_lists
  use libsupermesh_transform_elements
  use libsupermesh_elements
  use libsupermesh_construction
  use libsupermesh_tri_intersection_module, only : tri_buf_size
  use libsupermesh_tet_intersection_module, only : tet_buf_size

  type(vector_field) :: positionsA, positionsB
  type(ilist), dimension(:), allocatable :: map_BA
  real, dimension(:), allocatable :: detwei
  integer :: ele_A, ele_B, ele_C, n_trisC, n_tetsC, i
  real :: vol_B, vols_C
  logical :: fail
  type(inode), pointer :: llnode
  type(vector_field) :: intersection
  real, dimension(2, 3, tri_buf_size) ::  trisC_real
  real, dimension(3, 4, tet_buf_size) ::  tetsC_real
  type(mesh_type) :: new_mesh, intersection_mesh
  type(element_type) :: shape_

  positionsA = read_triangle_files("data/cube.1", quad_degree=6)
  positionsB = read_triangle_files("data/cube.2", quad_degree=6)

  allocate(map_BA(ele_count(positionsB)))
  allocate(detwei(ele_ngi(positionsA, 1)))

  map_BA = libsupermesh_advancing_front_intersection_finder( &
      & positionsB%val, reshape(positionsB%mesh%ndglno, (/ele_loc(positionsB, 1), ele_count(positionsB)/)), &
      & positionsA%val, reshape(positionsA%mesh%ndglno, (/ele_loc(positionsA, 1), ele_count(positionsA)/)) )
  call libsupermesh_intersector_set_dimension(positionsA%dim)
#ifdef HAVE_CGAL
  call libsupermesh_intersector_set_exactness(.true.)
#endif

  do ele_B=1,ele_count(positionsB)
    call transform_to_physical(positionsB, ele_B, detwei=detwei)
    vol_B = sum(detwei)

    llnode => map_BA(ele_B)%firstnode
    vols_C = 0.0
    do while(associated(llnode))
      ele_A = llnode%value
!      intersection = intersect_elements(positionsA, ele_A, ele_val(positionsB, ele_B), ele_shape(positionsB, ele_B))
      if ( positionsA%dim == 1 ) then
      ! 1D
        shape_ = ele_shape(positionsA, ele_A)
        intersection = libsupermesh_intersect_elements_old(ele_val(positionsA, ele_A), &
          ele_val(positionsB, ele_B), ele_loc(positionsA, ele_A), positionsA%dim, node_count(positionsA), &
          shape_%quadrature%vertices, shape_%quadrature%dim, shape_%quadrature%ngi, &
          shape_%quadrature%degree, shape_%loc, shape_%dim, shape_%degree)
      else if ( positionsA%dim == 2 ) then
        if ( ele_loc(positionsA, ele_A) == 3 ) then
        ! Triangles (2D)
          call libsupermesh_intersect_elements(ele_val(positionsA, ele_A), ele_val(positionsB, ele_B), &
                     n_trisC, trisC_real)
          call allocate(new_mesh, n_trisC * 3, n_trisC, ele_shape(positionsA, ele_A))
 
          if ( n_trisC > 0 ) then
            new_mesh%ndglno = (/ (i, i=1,3 * n_trisC) /)
            new_mesh%continuity = -1
          end if
 
          call allocate(intersection, positionsA%dim, new_mesh, "IntersectionCoordinates")
          if ( n_trisC > 0 ) then
            do i = 1, n_trisC
              call set(intersection, ele_nodes(intersection, i), trisC_real(:,:,i))
            end do
          end if
          call deallocate(new_mesh)
        else if ( ele_loc(positionsA, ele_A) == 4 ) then
        ! Quads (2D)
          shape_ = ele_shape(positionsA, ele_A)
          intersection = libsupermesh_intersect_elements_old(ele_val(positionsA, ele_A), &
            ele_val(positionsB, ele_B), ele_loc(positionsA, ele_A), positionsA%dim, node_count(positionsA), &
            shape_%quadrature%vertices, shape_%quadrature%dim, shape_%quadrature%ngi, &
            shape_%quadrature%degree, shape_%loc, shape_%dim, shape_%degree)
        end if
      else if ( positionsA%dim == 3 ) then
      ! Tets (3D)
        call libsupermesh_intersect_elements(ele_val(positionsA, ele_A), ele_val(positionsB, ele_B), &
                n_tetsC, tetsC_real )
        call allocate(intersection_mesh, n_tetsC * 4, n_tetsC, ele_shape(positionsB, ele_B))
        intersection_mesh%continuity = -1

        if ( n_tetsC > 0 ) then
          intersection_mesh%ndglno = (/ (i, i=1,4 * n_tetsC) /)
        end if
 
        call allocate(intersection, positionsA%dim, intersection_mesh, "IntersectionCoordinatesTet")
        if ( n_tetsC > 0 ) then
          do i = 1, n_tetsC
            call set(intersection, ele_nodes(intersection, i), tetsC_real(:,:,i))
          end do
        end if
        call deallocate(intersection_mesh)
      end if

      do ele_C=1,ele_count(intersection)
        call transform_to_physical(intersection, ele_C, detwei=detwei)
        vols_C = vols_C + sum(detwei)
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
