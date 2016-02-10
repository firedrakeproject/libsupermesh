subroutine test_tet_intersector

  use libsupermesh_construction
  use libsupermesh_fields
  use libsupermesh_read_triangle
  use libsupermesh_tri_intersection
  use libsupermesh_tet_intersection_module
  use libsupermesh_unittest_tools
  
  implicit none

  type(vector_field) :: positionsA, positionsB
  integer :: ele_A, ele_B, ele_C
  real :: vol_fort, vol_fort_public, &
      &    vol_F_intersect_elements, vol_G_intersect_elements
  logical :: fail
  
  type(tet_type) :: tet_A, tet_B
  type(tet_type), dimension(tet_buf_size) :: tetsC
  type(plane_type), dimension(4) :: planes_B
  integer :: n_tetsC
  real, dimension(3, 4, tet_buf_size) ::  tetsC_real
  integer :: i
  type(vector_field) :: intersection
  type(mesh_type) :: new_mesh

  integer, parameter :: dim = 3, loc = 4

  positionsA = read_triangle_files("data/plcC", dim)
  positionsB = read_triangle_files("data/plcD", dim)

  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)

      tet_A%v = ele_val(positionsA, ele_A)
      tet_B%v = ele_val(positionsB, ele_B)
      tet_A%colours = -1
      tet_B%colours = -1
      planes_B = get_planes(tet_B)

      ! B. Use the libSuperMesh internal triangle intersector (using derived types as input)
      call intersect_tets(tet_A, tet_B, tetsC, n_tetsC)
      vol_fort = 0.0
      do ele_C=1,n_tetsC
        vol_fort = vol_fort + tetrahedron_volume(tetsC(ele_C)%v)
      end do

      ! C. Use the libSuperMesh internal triangle intersector (using only reals as input)
      call intersect_tets(tet_A%v, tet_B%v, tetsC_real, n_tetsC)
      vol_fort_public = 0.0
      do ele_C=1,n_tetsC
        vol_fort_public = vol_fort_public + tetrahedron_volume(tetsC_real(:,:,ele_C))
      end do
      
      ! F. Use the new intersect_elements and do NOT create vector field
      call intersect_elements(tet_A%v, tet_B%v, n_tetsC, tetsC_real)
      vol_F_intersect_elements = 0.0
      do ele_C=1,n_tetsC
        vol_F_intersect_elements = vol_F_intersect_elements + tetrahedron_volume(tetsC_real(:,:,ele_C))
      end do

      ! G. Use the new intersect_elements and DO create vector field
      call intersect_elements(tet_A%v, tet_B%v, n_tetsC, tetsC_real)
      call allocate(new_mesh, dim, n_tetsC * loc, n_tetsC, loc)

      if ( n_tetsC > 0 ) then
        new_mesh%ndglno = (/ (i, i=1,loc * n_tetsC) /)
        new_mesh%continuity = -1
      end if

      call allocate(intersection, dim, new_mesh)
      if ( n_tetsC > 0 ) then
        do i = 1, n_tetsC
          call set(intersection, ele_nodes(intersection, i), tetsC_real(:,:,i))
        end do
      end if
      call deallocate(new_mesh)
      vol_G_intersect_elements = 0.0
      do ele_C=1,ele_count(intersection)
        vol_G_intersect_elements = vol_G_intersect_elements + tetrahedron_volume(ele_val(intersection, ele_C))
      end do
      call deallocate(intersection)

      fail = (vol_fort_public .fne. vol_fort) &
      & .or. (vol_F_intersect_elements .fne. vol_fort) &
      & .or. (vol_G_intersect_elements .fne. vol_fort)
      call report_test("[tet_intersector volumes]", fail, .false., "Should give the same volumes of intersection")
      if ( fail .eqv. .TRUE. ) then
        write (*,*) "[tet_intersector volumes]", &
          " vol_fort:", vol_fort,", vol_fort_public:",vol_fort_public,&
          ", vol_F_intersect_elements:",vol_F_intersect_elements, &
          ", vol_G_intersect_elements:",vol_G_intersect_elements,"."
      end if
      
    end do
  end do
  
  call deallocate(positionsA)
  call deallocate(positionsB)

end subroutine test_tet_intersector
