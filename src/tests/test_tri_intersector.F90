subroutine test_tri_intersector

  use libsupermesh_construction
  use libsupermesh_fields
  use libsupermesh_read_triangle
  use libsupermesh_tri_intersection_module
  use libsupermesh_unittest_tools
  
  implicit none

  type(vector_field) :: positionsA, positionsB
  integer :: ele_A, ele_B, ele_C
  real :: area_D_libwm, area_A_libwm_intersect, area_B_fort, area_C_fort_public, &
     &  area_E_intersect_elements, area_F_intersect_elements, area_G_intersect_elements
  logical :: fail
  
  type(tri_type) :: triA, triB
  type(tri_type), dimension(tri_buf_size) :: trisC
  integer :: ntests, n_trisC, i, nonods, totele, n_trisC_pre
  real, dimension(2, tri_buf_size) :: nodesC
  integer, dimension(3, tri_buf_size) :: ndglnoC
  real, dimension(2, 3, tri_buf_size) ::  trisC_real
  type(vector_field) :: intersection
  type(mesh_type) :: intersection_mesh, new_mesh
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
      call libsupermesh_intersect_tris(triA%v, triB%v, nodesC, ndglnoC, n_trisC)
      area_A_libwm_intersect = 0.0
      do ele_C=1,n_trisC
        area_A_libwm_intersect = area_A_libwm_intersect + triangle_area(nodesC(:, ndglnoC(:, ele_C)))
      end do
      n_trisC_pre = n_trisC

      ! B. Use the libSuperMesh internal triangle intersector (using derived types as input)
      call libsupermesh_intersect_tris(triA, triB, trisC, n_trisC)
      area_B_fort = 0.0
      do ele_C=1,n_trisC
        area_B_fort = area_B_fort + triangle_area(trisC(ele_C)%v)
      end do
      
      if (n_trisC_pre .ne. n_trisC) then
        write(*,*) "LibWM returned:",n_trisC_pre,", and tris returned:",n_trisC,"."
        write(*,*) "Input triangles"
        write(*,*) "A: ele_A:",ele_A,", triA:",ele_val(positionsA, ele_A),"."
        write(*,*) "B: ele_B:",ele_B,", triB:",ele_val(positionsB, ele_B),"."
        write(*,*) "Results:"
        do ele_C=1,n_trisC_pre
          write(*,*) "libWM:",nodesC(:, ndglnoC(:, ele_C)),"."
        end do
        do ele_C=1,n_trisC
          write(*,*) "tris :",trisC(ele_C)%v,"."
        end do
        call exit(1)
      end if

      ! C. Use the libSuperMesh internal triangle intersector (using only reals as input)
      call libsupermesh_intersect_tris(triA%v, triB%v, trisC_real, n_trisC)
      area_C_fort_public = 0.0
      do ele_C=1,n_trisC
        area_C_fort_public = area_C_fort_public + triangle_area(trisC_real(:,:,ele_C))
      end do

      ! D. Use libWM directly and create temporary vector field
      call libsupermesh_cintersector_set_input(ele_val(positionsA, ele_A), triB%v, positionsA%dim, ele_loc(positionsA, ele_A))
      call libsupermesh_cintersector_drive
      call libsupermesh_cintersector_query(nonods, totele)
      call allocate(intersection_mesh, nonods, totele, ele_shape(positionsA, ele_A), "TestIntersectionMesh")
      intersection_mesh%continuity = -1
      call allocate(intersection, positionsA%dim, intersection_mesh, "TestIntersectionCoordinates")
      if (nonods > 0) then
        call libsupermesh_cintersector_get_output(nonods, totele, positionsA%dim, ele_loc(positionsA, ele_A), nodes_tmp, intersection_mesh%ndglno)
        do i = 1, positionsA%dim
          intersection%val(i,:) = nodes_tmp((i - 1) * nonods + 1:i * nonods)
        end do
      end if
      area_D_libwm = 0.0
      do ele_C=1,totele
        area_D_libwm = area_D_libwm + abs(simplex_volume(intersection, ele_C))
      end do
      call deallocate(intersection_mesh)
      call deallocate(intersection)

      ! E. Use libWM directly and create temporary vector field      
      shape_lib = ele_shape(positionsA, ele_A)
      intersection = libsupermesh_intersect_elements_old(triA%v, &
        triB%v, ele_loc(positionsA, ele_A), positionsA%dim, node_count(positionsA), &
        shape_lib%quadrature%vertices, shape_lib%quadrature%dim, shape_lib%quadrature%ngi, &
        shape_lib%quadrature%degree, shape_lib%loc, shape_lib%dim, shape_lib%degree)
      area_E_intersect_elements = 0.0
      do ele_C=1,ele_count(intersection)
        area_E_intersect_elements = area_E_intersect_elements + abs(simplex_volume(intersection, ele_C))
      end do
      call deallocate(intersection)

      ! F. Use the new libsupermesh_intersect_elements and do NOT create vector field 
      call libsupermesh_intersect_elements(triA%v, triB%v, positionsA%dim, n_trisC, trisC_real=trisC_real)
      area_F_intersect_elements = 0.0
      do ele_C=1,n_trisC
        area_F_intersect_elements = area_F_intersect_elements + triangle_area(trisC_real(:,:,ele_C))
      end do

      ! G. Use the new libsupermesh_intersect_elements and DO create vector field 
      call libsupermesh_intersect_elements(triA%v, triB%v, positionsA%dim, n_trisC, trisC_real=trisC_real)
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
      area_G_intersect_elements = 0.0
      do ele_C=1,ele_count(intersection)
        area_G_intersect_elements = area_G_intersect_elements + abs(simplex_volume(intersection, ele_C))
      end do
      call deallocate(intersection)

      fail = (area_A_libwm_intersect .fne. area_B_fort) &
         .OR. (area_E_intersect_elements .fne. area_C_fort_public ) &
         .OR. (area_E_intersect_elements .fne. area_B_fort) &
         .OR. (area_A_libwm_intersect .fne. area_D_libwm) &
         .OR. (area_A_libwm_intersect .fne. area_C_fort_public ) &
         .OR. (area_F_intersect_elements .fne. area_A_libwm_intersect ) &
         .OR. (area_F_intersect_elements .fne. area_D_libwm ) &
         .OR. (area_G_intersect_elements .fne. area_B_fort) &
         .OR. (area_G_intersect_elements .fne. area_D_libwm )
      call report_test("[tri_intersector areas]",fail, .false., "Should give the same areas of intersection")
      if ( fail .eqv. .TRUE. ) then
        write (*,*) "[tri_intersector areas] area_A_libwm_intersect:",area_A_libwm_intersect, &
          ", area_B_fort:", area_B_fort,", area_C_fort_public:",area_C_fort_public,&
          ", area_D_libwm:",area_D_libwm, &
          ", area_E_intersect_elements:",area_E_intersect_elements, &
          ", area_F_intersect_elements:",area_F_intersect_elements, &
          ", area_G_intersect_elements:",area_G_intersect_elements,"."
      end if

    end do
  end do

  call libsupermesh_cintersection_finder_reset(ntests)
  call finalise_libsupermesh()
  
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
