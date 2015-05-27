subroutine test_tet_intersector

  use libsupermesh_construction
  use libsupermesh_fields
  use libsupermesh_read_triangle
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
  type(element_type) :: shape_lib
  
  ! I hope this is big enough ...
  real, dimension(1024) :: nodes_tmp

  positionsA = read_triangle_files("data/plcC",quad_degree=4)!, mdim=3)
  positionsB = read_triangle_files("data/plcD",quad_degree=4)!, mdim=3)

  call libsupermesh_intersector_set_dimension(3)
  call libsupermesh_intersector_set_exactness(.false.)

  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)

      tet_A%v = ele_val(positionsA, ele_A)
      tet_B%v = ele_val(positionsB, ele_B)
      tet_A%colours = -1
      tet_B%colours = -1
      planes_B = get_planes(tet_B)

      ! A. Use libWM without creating temporary vector fields.
      call libsupermesh_intersect_tets_libwm(tet_A%v, tet_B%v, nodesC, ndglnoC, n_tetsC)
      vol_libwm_intersect = 0.0
      do ele_C=1,n_tetsC
        vol_libwm_intersect = vol_libwm_intersect + abs(tetvol_test(nodesC(:, ndglnoC(:, ele_C))))
      end do

      ! B. Use the libSuperMesh internal triangle intersector (using derived types as input)
      call libsupermesh_intersect_tets(tet_A, tet_B, tetsC, n_tetsC)
      vol_fort = 0.0
      do ele_C=1,n_tetsC
        vol_fort = vol_fort + abs(tetvol_test(tetsC(ele_C)%v))
      end do

      ! C. Use the libSuperMesh internal triangle intersector (using only reals as input)
      call libsupermesh_intersect_tets(tet_A%v, tet_B%v, tetsC_real, n_tetsC)
      vol_fort_public = 0.0
      do ele_C=1,n_tetsC
        vol_fort_public = vol_fort_public + abs(tetvol_test(tetsC_real(:,:,ele_C)))
      end do

      ! D. Use libWM directly and create temporary vector field
      call libsupermesh_cintersector_set_input(ele_val(positionsA, ele_A), tet_B%v, positionsA%dim, ele_loc(positionsA, ele_A))
      call libsupermesh_cintersector_drive
      call libsupermesh_cintersector_query(nonods, totele)
      call allocate(intersection_mesh, nonods, totele, ele_shape(positionsA, ele_A), "IntersectionMesh")
      intersection_mesh%continuity = -1
      call allocate(intersection, positionsA%dim, intersection_mesh, "IntersectionCoordinates")
      if (nonods > 0) then
        call libsupermesh_cintersector_get_output(nonods, totele, positionsA%dim, ele_loc(positionsA, ele_A), nodes_tmp, intersection_mesh%ndglno)
        do i = 1, positionsA%dim
          intersection%val(i,:) = nodes_tmp((i - 1) * nonods + 1:i * nonods)
        end do
      end if
      vol_libwm = 0.0
      do ele_C=1,totele
        vol_libwm = vol_libwm + abs(simplex_volume(intersection, ele_C))
      end do
      call deallocate(intersection_mesh)
      call deallocate(intersection)

      ! E. Use libWM directly and create temporary vector field      
      shape_lib = ele_shape(positionsB, ele_B)
      intersect_elements_result = libsupermesh_intersect_elements_old(tet_B%v, &
        ele_val(positionsA, ele_A), &
        ele_loc(positionsB, ele_B), positionsB%dim, &
        node_count(positionsA), &
        shape_lib%quadrature%vertices, shape_lib%quadrature%dim, shape_lib%quadrature%ngi, &
        shape_lib%quadrature%degree, shape_lib%loc, shape_lib%dim, shape_lib%degree)
        
      vol_E_intersect_elements = 0.0
      do ele_C=1,ele_count(intersect_elements_result)
        vol_E_intersect_elements = vol_E_intersect_elements + abs(simplex_volume(intersect_elements_result, ele_C))
      end do
      call deallocate(intersect_elements_result)
      
      ! F. Use the new libsupermesh_intersect_elements and do NOT create vector field 
      call libsupermesh_intersect_elements(tet_A%v, tet_B%v, positionsA%dim, n_tetsC, tetsC_real=tetsC_real)
      vol_F_intersect_elements = 0.0
      do ele_C=1,n_tetsC
        vol_F_intersect_elements = vol_F_intersect_elements + abs(tetvol_test(tetsC_real(:,:,ele_C)))
      end do

      ! G. Use the new libsupermesh_intersect_elements and DO create vector field 
      call libsupermesh_intersect_elements(tet_A%v, tet_B%v, positionsA%dim, n_tetsC, tetsC_real=tetsC_real)
      call allocate(new_mesh, n_tetsC * 4, n_tetsC, ele_shape(positionsA, ele_A))

      if ( n_tetsC > 0 ) then
        new_mesh%ndglno = (/ (i, i=1,4 * n_tetsC) /)
        new_mesh%continuity = -1
      end if

      call allocate(intersection, positionsA%dim, new_mesh, "IntersectionCoordinates")
      if ( n_tetsC > 0 ) then
        do i = 1, n_tetsC
          call set(intersection, ele_nodes(intersection, i), tetsC_real(:,:,i))
        end do
      end if
      call deallocate(new_mesh)
      vol_G_intersect_elements = 0.0
      do ele_C=1,ele_count(intersection)
        vol_G_intersect_elements = vol_G_intersect_elements + abs(simplex_volume(intersection, ele_C))
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

contains

  function tetvol_test(positions) result(t)
    real, dimension(3, 4), intent(in) :: positions
    real :: t

    t = tetvol_test_old(positions(1, :), positions(2, :), positions(3, :))
    
  end function tetvol_test
  
  real function tetvol_test_old( x, y, z )

    real x(4), y(4), z(4)
    real vol, x12, x13, x14, y12, y13, y14, z12, z13, z14
    !
    x12 = x(2) - x(1)
    x13 = x(3) - x(1)
    x14 = x(4) - x(1)
    y12 = y(2) - y(1)
    y13 = y(3) - y(1)
    y14 = y(4) - y(1)
    z12 = z(2) - z(1)
    z13 = z(3) - z(1)
    z14 = z(4) - z(1)
    !
    vol = x12*( y13*z14 - y14*z13 )  &
         + x13*( y14*z12 - y12*z14 ) &
         + x14*( y12*z13 - y13*z12 )
    !
    tetvol_test_old = vol/6
    !
    return
  end function tetvol_test_old

end subroutine test_tet_intersector