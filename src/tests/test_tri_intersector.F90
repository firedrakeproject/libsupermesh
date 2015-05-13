subroutine test_tri_intersector

  use libsupermesh_read_triangle
  use libsupermesh_tri_intersection_module
  use libsupermesh_construction
  use libsupermesh_fields_base, only:ele_shape
  use libsupermesh_quadrature, only:make_quadrature
  use libsupermesh_shape_functions, only:make_element_shape
!  use libsupermesh_fields
  use libsupermesh_unittest_tools
  implicit none

  type(vector_field) :: positionsA, positionsB
  type(vector_field) :: libwm, fort
  integer :: ele_A, ele_B, ele_C
  real :: area_libwm, area_fort
  logical :: fail
  integer :: stat
  type(tri_type) :: tri_A, tri_B
  type(line_type), dimension(3) :: lines_B
  
!  real, dimension(:,:), allocatable :: positions_B_lib_val
  real, dimension(:,:), allocatable :: positions_B_lib_val
  type(quadrature_type) :: quad_lib
  type(element_type) :: shape_lib
  integer :: dimB, n_count

  positionsA = read_triangle_files("data/plcA", quad_degree=3, mdim=2)
  positionsB = read_triangle_files("data/plcB", quad_degree=3, mdim=2)

  call libsupermesh_intersector_set_dimension(2)
  call libsupermesh_intersector_set_exactness(.false.)

  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)
      
      dimB = positionsB%dim
      n_count = 0
      select case(positionsB%field_type)
      case(FIELD_TYPE_NORMAL)
        n_count = node_count(positionsB%mesh)
!        allocate(positions_B_lib_val(dimB,n_count))
      case(FIELD_TYPE_CONSTANT)
!        allocate(positions_B_lib_val(dimB,1))
      case(FIELD_TYPE_DEFERRED)
!        allocate(positions_B_lib_val(0,0))
      end select
       allocate(positions_B_lib_val(dimB, dimB+1))
      positions_B_lib_val = ele_val(positionsB, ele_B)
      shape_lib = ele_shape(positionsB, ele_B)

      libwm = libsupermesh_intersect_elements(positions_B_lib_val, &
        ele_val(positionsA, ele_A), &
        ele_loc(positionsB, ele_B), dimB, &
        node_count(positionsA), &
        shape_lib%quadrature%vertices, shape_lib%quadrature%dim, shape_lib%quadrature%ngi, &
        shape_lib%quadrature%degree, shape_lib%loc, shape_lib%dim, shape_lib%degree)
      deallocate(positions_B_lib_val)
      tri_A%v = ele_val(positionsA, ele_A)
      tri_B%v = ele_val(positionsB, ele_B)
      lines_B = get_lines(tri_B)
      call libsupermesh_intersect_tris(tri_A, lines_B, shape=ele_shape(positionsB, 1), stat=stat, output=fort)
      fail = (ele_count(libwm) /= ele_count(fort))
!      call report_test("[tet_intersector counts]", fail, .false., "Should give the same number of elements.")

      area_libwm = 0.0
      do ele_C=1,ele_count(libwm)
        area_libwm = area_libwm + abs(simplex_volume(libwm, ele_C))
!        write(*,*) "libsupermesh_intersect_elements: ",ele_val(libwm, ele_C),"."
      end do
      area_fort = 0.0
      if (stat == 0) then
        do ele_C=1,ele_count(fort)
          area_fort = area_fort + abs(simplex_volume(fort, ele_C))
        end do
      end if
      fail = (area_libwm .fne. area_fort)
      call report_test("[tri_intersector areas]", fail, .false., "Should give the same areas of intersection")
      if ( fail .eqv. .TRUE. ) then
        write (*,*) "[tri_intersector areas] area_libwm:",area_libwm,", area_fort:",area_fort,"."
      end if

      call deallocate(libwm)
      if (stat == 0) then
        call deallocate(fort)
      end if
    end do
  end do
  call deallocate(positionsA)
  call deallocate(positionsB)

end subroutine test_tri_intersector