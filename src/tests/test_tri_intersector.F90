subroutine test_tri_intersector

  use libsupermesh_construction
  use libsupermesh_fields_allocates
  use libsupermesh_fields_base
  use libsupermesh_read_triangle
  use libsupermesh_tri_intersection_module
  use libsupermesh_unittest_tools
  
  implicit none

  type(vector_field) :: positionsA, positionsB
  integer :: ele_A, ele_B, ele_C
  real :: area_libwm, area_libwm_intersect, area_fort, area_fort_public, area_intersect_elements
  logical :: fail
  
  type(tri_type) :: triA, triB
  type(tri_type), dimension(8) :: trisC
  integer :: ntests, n_trisC, i, nonods, totele
  real, dimension(2, 8) :: nodesC
  integer, dimension(3, 8) :: ndglnoC
  real, dimension(2, 3, 8) ::  trisC_real
  type(vector_field) :: intersection
  type(mesh_type) :: intersection_mesh
  type(element_type) :: shape_lib

  ! I hope this is big enough ...
  real, dimension(1024) :: nodes_tmp
  
  positionsA = read_triangle_files("data/plcA", quad_degree=0, mdim=2)
  positionsB = read_triangle_files("data/plcB", quad_degree=0, mdim=2)

  call libsupermesh_intersector_set_dimension(2)
  call libsupermesh_intersector_set_exactness(.false.)

  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)
      
      triA%v = ele_val(positionsA, ele_A)
      triB%v = ele_val(positionsB, ele_B)

      ! A. Use libWM without creating temporary vector fields.
      call libsupermesh_intersect_tris_libwm(triA%v, triB%v, nodesC, ndglnoC, n_trisC)
      area_libwm_intersect = 0.0
      do ele_C=1,n_trisC
        area_libwm_intersect = area_libwm_intersect + triangle_area(nodesC(:, ndglnoC(:, ele_C)))
      end do

      ! B. Use the libSuperMesh internal triangle intersector (using derived types as input)
      call libsupermesh_intersect_tris_dt(triA, triB, trisC, n_trisC)
      area_fort = 0.0
      do ele_C=1,n_trisC
        area_fort = area_fort + triangle_area(trisC(ele_C)%v)
      end do

      ! C. Use the libSuperMesh internal triangle intersector (using only reals as input)
      call libsupermesh_intersect_tris_dt_public(triA%v, triB%v, trisC_real, n_trisC)
      area_fort_public = 0.0
      do ele_C=1,n_trisC
        area_fort_public = area_fort_public + triangle_area(trisC_real(:,:,ele_C))
      end do

      ! D. Use libWM directly and create temporary vector field
      call libsupermesh_cintersector_set_input(ele_val(positionsA, ele_A), triB%v, positionsA%dim, ele_loc(positionsA, ele_A))
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
      area_libwm = 0.0
      do ele_C=1,totele
        area_libwm = area_libwm + abs(simplex_volume(intersection, ele_C))
      end do
      call deallocate(intersection_mesh)
      call deallocate(intersection)

      ! E. Use libWM directly and create temporary vector field      
      shape_lib = ele_shape(positionsA, ele_A)
      intersection = libsupermesh_intersect_elements(triA%v, &
        triB%v, ele_loc(positionsA, ele_A), positionsA%dim, node_count(positionsA), &
        shape_lib%quadrature%vertices, shape_lib%quadrature%dim, shape_lib%quadrature%ngi, &
        shape_lib%quadrature%degree, shape_lib%loc, shape_lib%dim, shape_lib%degree)
      area_intersect_elements = 0.0
      do ele_C=1,ele_count(intersection)
        area_intersect_elements = area_intersect_elements + abs(simplex_volume(intersection, ele_C))
        area_intersect_elements = area_intersect_elements
      end do
      call deallocate(intersection)
      
      fail = (area_libwm_intersect .fne. area_fort) &
         .OR. (area_intersect_elements .fne. area_fort_public ) &
         .OR. (area_intersect_elements .fne. area_fort) &
         .OR. (area_libwm_intersect .fne. area_libwm) &
         .OR. (area_libwm_intersect .fne. area_fort_public )
      call report_test("[tri_intersector areas]",fail, .false., "Should give the same areas of intersection")
      if ( fail .eqv. .TRUE. ) then
        write (*,*) "[tri_intersector areas] area_libwm_intersect:",area_libwm_intersect, &
          ", area_fort:", area_fort,", area_fort_public:",area_fort_public,&
          ", area_intersect_elements:",area_intersect_elements,"."
      end if

    end do
  end do

  call libsupermesh_cintersection_finder_reset(ntests)
  
  call deallocate(positionsA)
  call deallocate(positionsB)

contains

  pure function triangle_area(tri) result(area)
    real, dimension(2, 3), intent(in) :: tri

    real :: area
    real, dimension(2) :: u, v

    u = tri(:, 3) - tri(:, 1)
    v = tri(:, 2) - tri(:, 1)

    area = 0.5 * abs(u(2) * v(1) - u(1) * v(2))

  end function triangle_area

end subroutine test_tri_intersector