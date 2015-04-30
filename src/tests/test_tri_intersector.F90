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
  real :: vol_libwm, vol_fort
  logical :: fail
  integer :: stat
  type(tri_type) :: tri_A, tri_B
  type(line_type), dimension(3) :: lines_B
  
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
        allocate(positions_B_lib_val(dimB,n_count))
      case(FIELD_TYPE_CONSTANT)
        allocate(positions_B_lib_val(dimB,1))
      case(FIELD_TYPE_DEFERRED)
        allocate(positions_B_lib_val(0,0))
      end select
      positions_B_lib_val = positionsB%val
      
      write (*,*) "test_tri_intersector: Begin."
        
!      quad_lib = make_quadrature(vertices = positionsB%mesh%shape%quadrature%vertices, dim = positionsB%mesh%shape%quadrature%dim, ngi = positionsB%mesh%shape%quadrature%ngi, degree = positionsB%mesh%shape%degree * 2)
!      shape_lib = make_element_shape(vertices = positionsB%mesh%shape%loc, dim = positionsB%mesh%shape%dim, degree = positionsB%mesh%shape%degree, quad = quad_lib)
!      call deallocate(quad_lib)
      
      write (*,*) "test_tri_intersector: positionsA%field_type:",positionsA%field_type,&
        &", positionsA%dim:",positionsA%dim,", positionsA%shape%loc:",positionsA%mesh%shape%loc,"."
      
!      libwm = libsupermesh_intersect_elements(positionsB, ele_B, ele_val(positionsA, ele_A), ele_shape(positionsB, 1))
      libwm = libsupermesh_intersect_elements(positions_B_lib_val, ele_count(positionsB), &
        positionsB%mesh%shape%quadrature%vertices, positionsB%mesh%shape%quadrature%dim, ele_B, &
        ele_val(positionsA, ele_A), ele_shape(positionsB, 1), ele_loc(positionsB, ele_B), dimB, &
        node_count(positionsB), positionsB%mesh%shape%loc, positionsB%field_type, positionsB%mesh%ndglno)
!      call deallocate(shape_lib)
      deallocate(positions_B_lib_val)
      write (*,*) "test_tri_intersector:1, ele_A:",ele_A,"."
      tri_A%v = ele_val(positionsA, ele_A)
      write (*,*) "test_tri_intersector:2, ele_B:",ele_B,"."
      tri_B%v = ele_val(positionsB, ele_B)
      write (*,*) "test_tri_intersector:3"
      lines_B = get_lines(tri_B)
      write (*,*) "test_tri_intersector:4"
      call libsupermesh_intersect_tris(tri_A, lines_B, shape=ele_shape(positionsB, 1), stat=stat, output=fort)
      write (*,*) "test_tri_intersector:5"
      fail = (ele_count(libwm) /= ele_count(fort))
      write (*,*) "number of elements (libwm:",ele_count(libwm),", fort:",ele_count(fort),"."
!      call report_test("[tet_intersector counts]", fail, .false., "Should give the same number of elements.")

      vol_libwm = 0.0
      do ele_C=1,ele_count(libwm)
        vol_libwm = vol_libwm + abs(simplex_volume(libwm, ele_C))
      end do
      vol_fort = 0.0
      if (stat == 0) then
        do ele_C=1,ele_count(fort)
          vol_fort = vol_fort + abs(simplex_volume(fort, ele_C))
        end do
      end if
      write (*,*) "test_tri_intersector:6"
      fail = (vol_libwm .fne. vol_fort)
      write (*,*) "Should give the same area (libwm:",vol_libwm,", fort:",vol_fort,")."
      call report_test("[tri_intersector areas]", fail, .false., "Should give the same areas of intersection")

      call deallocate(libwm)
      write (*,*) "test_tri_intersector:7"
      if (stat == 0) then
        call deallocate(fort)
      end if
      write (*,*) "test_tri_intersector:8"
    end do
  end do
  write (*,*) "test_tri_intersector:9"
  call deallocate(positionsA)
  write (*,*) "test_tri_intersector:10"
  call deallocate(positionsB)
  write (*,*) "test_tri_intersector:11"

end subroutine test_tri_intersector