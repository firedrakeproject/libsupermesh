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
  type(plane_type), dimension(4) :: planes_B
  
  real, dimension(:,:), allocatable :: positions_B_lib_val
  type(quadrature_type) :: quad_lib
  type(element_type) :: shape_lib
  integer :: dimB, n_count

  positionsA = read_triangle_files("data/plcA", quad_degree=4)
  positionsB = read_triangle_files("data/plcB", quad_degree=4)

  call libsupermesh_intersector_set_dimension(3)
  call libsupermesh_intersector_set_exactness(.false.)

  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)
      
      dimB = positionsB%dim
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
      
      quad_lib = make_quadrature(vertices = positionsB%mesh%shape%quadrature%vertices, dim = positionsB%mesh%shape%quadrature%dim, ngi = positionsB%mesh%shape%quadrature%ngi, degree = positionsB%mesh%shape%degree * 2)
      shape_lib = make_element_shape(vertices = positionsB%mesh%shape%loc, dim = positionsB%mesh%shape%dim, degree = positionsB%mesh%shape%degree, quad = quad_lib)
      call deallocate(quad_lib)
      
      write (*,*) "test_tri_intersector: positionsA%field_type:",positionsA%field_type,&
        &", positionsA%dim:",positionsA%dim,", positionsA%shape%loc:",positionsA%mesh%shape%loc,"."
      
!      libwm = libsupermesh_intersect_elements(positionsB, ele_B, ele_val(positionsA, ele_A), ele_shape(positionsB, 1))
      libwm = libsupermesh_intersect_elements(positions_B_lib_val, ele_count(positionsB), &
        positionsB%mesh%shape%quadrature%vertices, positionsB%mesh%shape%quadrature%dim, ele_B, &
        ele_val(positionsA, ele_A), shape_lib, ele_loc(positionsB, ele_B), dimB, node_count(positionsB), &
        positionsB%mesh%shape%loc, positionsB%field_type, positionsB%mesh%ndglno)
      call deallocate(shape_lib)
      write (*,*) "test_tri_intersector:1"
      tri_A%v = ele_val(positionsA, ele_A)
      write (*,*) "test_tri_intersector:2"
      tri_B%v = ele_val(positionsB, ele_B)
      write (*,*) "test_tri_intersector:3"
      planes_B = get_planes(tri_B)
      call libsupermesh_intersect_tris(tri_A, planes_B, shape=ele_shape(positionsB, 1), stat=stat, output=fort)

      fail = (ele_count(libwm) /= ele_count(fort))
!      call report_test("[tet_intersector counts]", fail, .false., "Should give the same number of elements")

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

!      fail = (vol_libwm .fne. vol_fort)
!      call report_test("[tet_intersector volumes]", fail, .false., "Should give the same volumes of intersection")

      call deallocate(libwm)
      if (stat == 0) then
        call deallocate(fort)
      end if
    end do
  end do
  call deallocate(positionsA)
  call deallocate(positionsB)

end subroutine test_tri_intersector