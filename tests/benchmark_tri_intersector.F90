#define BUF_SIZE_A 50
#define BUF_SIZE_B 2000
subroutine benchmark_tri_intersector

  use libsupermesh_read_triangle
  use libsupermesh_tri_intersection_module
  use libsupermesh_construction
  use libsupermesh_fields_base, only:ele_shape
  use libsupermesh_quadrature, only:make_quadrature
  use libsupermesh_shape_functions, only:make_element_shape
  use libsupermesh_unittest_tools
  implicit none
  
  include "mpif.h"

  type(vector_field) :: positionsAsmall, positionsBsmall, positionsA, positionsB
  type(vector_field) :: libwm, fort
  integer :: ele_A, ele_B, ele_C
  real, dimension(BUF_SIZE_A * BUF_SIZE_B) :: area_libwm, area_fort
  real :: local_area_libwm, local_area_fort
  logical :: fail, totalFail = .FALSE.
  integer :: stat
  type(tri_type) :: tri_A, tri_B
  type(line_type), dimension(3) :: lines_B
  
  real, dimension(:,:), allocatable :: positions_B_lib_val
  type(quadrature_type) :: quad_lib
  type(element_type) :: shape_lib
  integer :: dimB, n_count, ele, i, index, elementCount, loc, nodeCount
  
  type(mesh_type) :: intersection_mesh_A, intersection_mesh_B
  type(tri_type), dimension(BUF_SIZE_A) :: local_tri_array
  type(element_type) :: shape
  REAL, dimension(2,3) :: num
  double precision :: t1,t2,dtWM,dtLibSuperMesh
  real, dimension(:, :), allocatable :: posB
  type(element_type) :: elementShape

  positionsAsmall = read_triangle_files("data/plcA", quad_degree=3, mdim=2)
  positionsBsmall = read_triangle_files("data/plcB", quad_degree=3, mdim=2)

  call libsupermesh_intersector_set_dimension(2)
  call libsupermesh_intersector_set_exactness(.false.)
  
  CALL RANDOM_SEED ()
  
  call allocate(intersection_mesh_A, BUF_SIZE_A * 3, BUF_SIZE_A, ele_shape(positionsAsmall, 1), name="IntersectionMeshA")
  intersection_mesh_A%ndglno = (/ (i, i=1,BUF_SIZE_A*3) /)
  intersection_mesh_A%continuity = -1
      
  intersection_mesh_A%nodes = BUF_SIZE_A*3
  intersection_mesh_A%elements = BUF_SIZE_A
  call allocate(positionsA, 2, intersection_mesh_A, "IntersectionCoordinatesA")
  
  do ele=1,BUF_SIZE_A
    CALL RANDOM_NUMBER(num)
    num = num * 10
    call set(positionsA, ele_nodes(positionsA, ele), num)
  end do
  
  call allocate(intersection_mesh_B, BUF_SIZE_B * 3, BUF_SIZE_B, ele_shape(positionsBsmall, 1), name="IntersectionMeshB")
  intersection_mesh_B%ndglno = (/ (i, i=1,BUF_SIZE_B*3) /)
  intersection_mesh_B%continuity = -1
      
  intersection_mesh_B%nodes = BUF_SIZE_B*3
  intersection_mesh_B%elements = BUF_SIZE_B
  call allocate(positionsB, 2, intersection_mesh_B, "IntersectionCoordinatesB")
  
  do ele=1,BUF_SIZE_B
    CALL RANDOM_NUMBER(num)
    num = num * 10
    call set(positionsB, ele_nodes(positionsB, ele), num)
  end do
  
  dtWM = 0.0
  dtLibSuperMesh = 0.0
  area_libwm = 0.0
  area_fort = 0.0

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
      
      elementCount = ele_count(positionsB)
      posB = ele_val(positionsA, ele_A)
      elementShape = ele_shape(positionsB, 1)
      loc = ele_loc(positionsB, ele_B)
      nodeCount = node_count(positionsB)

      t1 = MPI_Wtime();
      libwm = libsupermesh_intersect_elements(positions_B_lib_val, elementCount, &
        positionsB%mesh%shape%quadrature%vertices, positionsB%mesh%shape%quadrature%dim, ele_B, &
        posB, elementShape, loc, dimB, &
        nodeCount, positionsB%mesh%shape%loc, positionsB%field_type, positionsB%mesh%ndglno)
      t2 = MPI_Wtime();
      dtWM = dtWM + ( t2 -t1 )
      
      deallocate(positions_B_lib_val)

      index = (ele_A-1)*BUF_SIZE_B + ele_B
      do ele_C=1,ele_count(libwm)
        area_libwm(index)  = area_libwm(index) + abs(simplex_volume(libwm, ele_C))
      end do
!      write (*,*) "benchmark_tri_intersector: area_libwm(",index,"):",area_libwm(index),"."

      call deallocate(libwm)
    end do
  end do
  write (*,*) "benchmark_tri_intersector: "
  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)
      
      tri_A%v = ele_val(positionsA, ele_A)
      tri_B%v = ele_val(positionsB, ele_B)
      t1 = MPI_Wtime();
      lines_B = get_lines(tri_B)
      call libsupermesh_intersect_tris(tri_A, lines_B, shape=ele_shape(positionsB, 1), stat=stat, output=fort)
      t2 = MPI_Wtime();
      dtLibSuperMesh = dtLibSuperMesh + ( t2 -t1 )

      index = (ele_A-1)*BUF_SIZE_B + ele_B
!      if ( index .eq. 31747 ) then
!        write (*,*) "benchmark_tri_intersector: triA:",tri_A,", triB:",tri_B,"."
!      end if
      
      if (stat == 0) then
        do ele_C=1,ele_count(fort)
          area_fort(index) = area_fort(index) + abs(simplex_volume(fort, ele_C))
        end do
      end if
!      write (*,*) "benchmark_tri_intersector: area_fort(",index,"):",area_fort(index),"."

      if (stat == 0) then
        call deallocate(fort)
      end if
    end do
  end do
  
  write (*,*) "benchmark_tri_intersector: dtWM:",dtWM,", dtLibSuperMesh:",dtLibSuperMesh,"."
  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)
      index = (ele_A-1)*BUF_SIZE_B + ele_B
      fail = (area_libwm(index) .fne. area_fort(index))
      
      if ( fail ) then
        write (*,*) "benchmark_tri_intersector: index:",index," area_libwm:",area_libwm(index),", area_fort:",area_fort(index),"."
        totalFail = .TRUE.
        call report_test("[benchmark_tri_intersector areas]", fail, .false., "Should give the same areas of intersection")
      end if
    end do
  end do
  
  call report_test("[benchmark_tri_intersector areas]", totalFail, .false., "Should give the same areas of intersection")
  
  call deallocate(positionsAsmall)
  call deallocate(positionsBsmall)
  call deallocate(intersection_mesh_A)
  call deallocate(intersection_mesh_B)
  call deallocate(positionsA)
  call deallocate(positionsB)

end subroutine benchmark_tri_intersector