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
  real :: area_libwm, area_fort
  logical :: fail
  
  type(tri_type) :: triA, triB
  type(tri_type), dimension(8) :: trisC
  integer :: ntests, n_trisC
  real, dimension(2, 8) :: nodesC
  integer, dimension(3, 8) :: ndglnoC

  positionsA = read_triangle_files("data/plcA", quad_degree=0, mdim=2)
  positionsB = read_triangle_files("data/plcB", quad_degree=0, mdim=2)

  call libsupermesh_intersector_set_dimension(2)
  call libsupermesh_intersector_set_exactness(.false.)

  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)
      
      triA%v = ele_val(positionsA, ele_A)
      triB%v = ele_val(positionsB, ele_B)

      call libsupermesh_intersect_tris_libwm(triA%v, triB%v, nodesC, ndglnoC, n_trisC)
      area_libwm = 0.0
      do ele_C=1,n_trisC
        area_libwm = area_libwm + triangle_area(nodesC(:, ndglnoC(:, ele_C)))
      end do
      
      call libsupermesh_intersect_tris_dt(triA, triB, trisC, n_trisC)
      area_fort = 0.0
      do ele_C=1,n_trisC
        area_fort = area_fort + triangle_area(trisC(ele_C)%v)
      end do
      fail = (area_libwm .fne. area_fort)
      call report_test("[tri_intersector areas]", fail, .false., "Should give the same areas of intersection")
      if ( fail .eqv. .TRUE. ) then
        write (*,*) "[tri_intersector areas] area_libwm:",area_libwm,", area_fort:",area_fort,"."
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