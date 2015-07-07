subroutine test_tet_intersector

  use libsupermesh_construction
  use libsupermesh_fields_dummy
  use libsupermesh_read_triangle_2
  use libsupermesh_tri_intersection_module
  use libsupermesh_tet_intersection_module
  use libsupermesh_unittest_tools
  
  implicit none

  type(vector_field) :: positionsA, positionsB
  integer :: ele_A, ele_B, ele_C
  real :: vol_libwm, vol_libwm_intersect, vol_fort, vol_fort_public, vol_E_intersect_elements, &
      &    vol_F_intersect_elements, vol_G_intersect_elements
  logical :: fail
  
  type(tet_type) :: tet_A, tet_B
  type(tet_type), dimension(tet_buf_size) :: tetsC
  type(plane_type), dimension(4) :: planes_B
  integer :: ntests, n_tetsC
  real, dimension(3, tet_buf_size) :: nodesC
  integer, dimension(4, tet_buf_size) :: ndglnoC
  real, dimension(3, 4, tet_buf_size) ::  tetsC_real
  integer :: i, nonods, totele
  type(vector_field) :: intersection, intersect_elements_result
  type(mesh_type) :: intersection_mesh, new_mesh
  
  ! I hope this is big enough ...
  real, dimension(1024) :: nodes_tmp

  integer, parameter :: dim = 3, loc = 4

  positionsA = read_triangle_files("data/plcC", dim)
  positionsB = read_triangle_files("data/plcD", dim)

  call libsupermesh_cintersector_set_dimension(dim)

  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)

      tet_A%v = ele_val(positionsA, ele_A)
      tet_B%v = ele_val(positionsB, ele_B)
      tet_A%colours = -1
      tet_B%colours = -1
      planes_B = get_planes(tet_B)

      ! A. Use libWM without creating temporary vector fields.
      call intersect_tets_libwm(tet_A%v, tet_B%v, nodesC, ndglnoC, n_tetsC)
      vol_libwm_intersect = 0.0
      do ele_C=1,n_tetsC
        vol_libwm_intersect = vol_libwm_intersect + tetrahedron_volume(nodesC(:, ndglnoC(:, ele_C)))
      end do

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

      ! D. Use libWM directly and create temporary vector field
      call libsupermesh_cintersector_set_input(ele_val(positionsA, ele_A), tet_B%v, dim, loc)
      call libsupermesh_cintersector_drive
      call libsupermesh_cintersector_query(nonods, totele)
      call allocate(intersection_mesh, dim, nonods, totele, loc)
      intersection_mesh%continuity = -1
      call allocate(intersection, dim, intersection_mesh)
      if (nonods > 0) then
        call libsupermesh_cintersector_get_output(nonods, totele, dim, loc, nodes_tmp, intersection_mesh%ndglno)
        do i = 1, dim
          intersection%val(i,:) = nodes_tmp((i - 1) * nonods + 1:i * nonods)
        end do
      end if
      vol_libwm = 0.0
      do ele_C=1,totele
        vol_libwm = vol_libwm + tetrahedron_volume(ele_val(intersection, ele_C))
      end do
      call deallocate(intersection_mesh)
      call deallocate(intersection)

      ! E. Use libWM directly and create temporary vector field
      intersect_elements_result = intersect_elements_old(tet_B%v, &
        ele_val(positionsA, ele_A), loc, dim)
      vol_E_intersect_elements = 0.0
      do ele_C=1,ele_count(intersect_elements_result)
        vol_E_intersect_elements = vol_E_intersect_elements + tetrahedron_volume(ele_val(intersect_elements_result, ele_C))
      end do
      call deallocate(intersect_elements_result)
      
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

      fail = (vol_libwm_intersect .fne. vol_fort) &
         .OR. (vol_E_intersect_elements .fne. vol_fort_public ) &
         .OR. (vol_E_intersect_elements .fne. vol_fort) &
         .OR. (vol_libwm_intersect .fne. vol_libwm) &
         .OR. (vol_libwm_intersect .fne. vol_fort_public ) &
         .OR. (vol_F_intersect_elements .fne. vol_libwm_intersect ) &
         .OR. (vol_F_intersect_elements .fne. vol_libwm ) &
         .OR. (vol_G_intersect_elements .fne. vol_fort) &
         .OR. (vol_G_intersect_elements .fne. vol_libwm )
      call report_test("[tet_intersector volumes]", fail, .false., "Should give the same volumes of intersection")
      if ( fail .eqv. .TRUE. ) then
        write (*,*) "[tet_intersector volumes] vol_libwm:",vol_libwm, &
          ", vol_fort:", vol_fort,", vol_fort_public:",vol_fort_public,&
          ", vol_E_intersect_elements:",vol_E_intersect_elements, &
          ", vol_F_intersect_elements:",vol_F_intersect_elements, &
          ", vol_G_intersect_elements:",vol_G_intersect_elements,"."
      end if
      
!    write(*,*) "test_tet_intersector: vol_libwm:",vol_libwm,", vol_fort:",vol_fort,", vol_fort_public:",vol_fort_public,", vol_intersect_elements:",vol_intersect_elements,"."
    end do
  end do

  call libsupermesh_cintersection_finder_reset(ntests)
  
  call deallocate(positionsA)
  call deallocate(positionsB)

end subroutine test_tet_intersector
