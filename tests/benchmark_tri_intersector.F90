#define BUF_SIZE_A 5
#define BUF_SIZE_B 30
subroutine benchmark_tri_intersector

  use libsupermesh_construction
  use libsupermesh_fields_allocates
  use libsupermesh_fields_base
  use libsupermesh_read_triangle
  use libsupermesh_tri_intersection_module
  use libsupermesh_unittest_tools
  
  implicit none
  
  include "mpif.h"

  type(vector_field) :: positionsA, positionsB
  integer :: ele_A, ele_B, ele_C
  real, dimension(BUF_SIZE_A * BUF_SIZE_B) :: area_libwm, area_fort, area_intersect, area_intersect_libwm
  logical :: fail1, fail2, fail3, fail4, totalFail = .FALSE.
  
  type(tri_type) :: triA, triB
  type(tri_type), dimension(8) :: trisC
  integer :: ntests, n_trisC
  real, dimension(2, 8) :: nodesC
  integer, dimension(3, 8) :: ndglnoC
  
  type(vector_field) :: libsupermesh_intersect_elements_result, libsupermesh_tri_intersector_result, libwm_result
  real :: local_area_libwm, local_area_fort
  double precision :: t1,t2,dtLibWM,dtLibSuperMeshTriIntersector,dtLibSuperMeshIntersectElements,dtLibSuperMeshIntersectLibWM
  real, dimension(:, :), allocatable :: posB, posA
  type(mesh_type) :: intersection_mesh_A, intersection_mesh_B, intersection_meshLibWM
  REAL, dimension(2) :: num

  integer :: index, ele, dimB, i, locB, loc, elementCount, n_count, nodeCount
  type(element_type) :: elementShape
  
  real, dimension(:,:), allocatable :: positions_B_lib_val
  
  integer :: nonods, totele
  
  type(mesh_type) :: mesh_lib
  
  ! I hope this is big enough ...
  real, dimension(1024) :: nodes_tmp

  call libsupermesh_intersector_set_dimension(2)
  call libsupermesh_intersector_set_exactness(.false.)
  
  open (unit = 20, file = "plcC_temp.node")
  open (unit = 21, file = "plcC_temp.ele")
  open (unit = 30, file = "plcD_temp.node")
  open (unit = 31, file = "plcD_temp.ele")
  
  CALL RANDOM_SEED ()
  
  write (20,*) "",BUF_SIZE_A*3," 2 0 0"
  write (21,*) "",BUF_SIZE_A," 3 0"
  do ele=1,BUF_SIZE_A*3
    CALL RANDOM_NUMBER(num)
    num = num * 10
    write (20,*) " ",ele," ",num," 0"
  end do
  do ele=1,BUF_SIZE_A
    write (21,*) " ",ele," ",(ele*3)-2," ", (ele*3)-1," ",(ele*3)
  end do
  
  write (30,*) "",BUF_SIZE_B*3," 2 0 0"
  write (31,*) "",BUF_SIZE_B," 3 0"
  do ele=1,BUF_SIZE_B*3
    CALL RANDOM_NUMBER(num)
    num = num * 10
    write (30,*) " ",ele," ",num," 0"
  end do
  do ele=1,BUF_SIZE_B
    write (31,*) " ",ele," ",(ele*3)-2," ", (ele*3)-1," ",(ele*3)
  end do
  
  close(20)
  close(21)
  close(30)
  close(31)
  
  positionsA = read_triangle_files("plcC_temp", quad_degree=0, mdim=2)
  positionsB = read_triangle_files("plcD_temp", quad_degree=0, mdim=2)
  
  dtLibWM = 0.0
  dtLibSuperMeshTriIntersector = 0.0
  dtLibSuperMeshIntersectElements = 0.0
  dtLibSuperMeshIntersectLibWM = 0.0
  area_libwm = 0.0
  area_fort = 0.0
  area_intersect = 0.0
  area_intersect_libwm = 0.0

  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)
       dimB = positionsB%dim
       n_count = 0
       select case(positionsB%field_type)
       case(FIELD_TYPE_NORMAL)
!         n_count = node_count(positionsB%mesh)
         n_count = positionsB%mesh%shape%quadrature%vertices
         allocate(positions_B_lib_val(dimB,n_count))
       case(FIELD_TYPE_CONSTANT)
         allocate(positions_B_lib_val(dimB,1))
       case(FIELD_TYPE_DEFERRED)
         allocate(positions_B_lib_val(0,0))
       end select
       positions_B_lib_val = ele_val(positionsB, ele_B)

       elementCount = ele_count(positionsB)
       allocate(posB(positionsA%dim, ele_loc(positionsA, ele_A)))
       posB = ele_val(positionsA, ele_A)
       elementShape = ele_shape(positionsB, ele_B)
       loc = ele_loc(positionsB, ele_B)
       nodeCount = node_count(positionsB)

       t1 = MPI_Wtime();
       libsupermesh_intersect_elements_result = libsupermesh_intersect_elements(positions_B_lib_val, &
         posB, loc, dimB, nodeCount, &
         elementShape%quadrature%vertices, elementShape%quadrature%dim, elementShape%quadrature%ngi, &
         elementShape%quadrature%degree, elementShape%loc, elementShape%dim, elementShape%degree)
       t2 = MPI_Wtime();
       dtLibSuperMeshIntersectElements = dtLibSuperMeshIntersectElements + ( t2 -t1 )

       deallocate(positions_B_lib_val)
       deallocate(posB)

       index = (ele_A-1)*BUF_SIZE_B + ele_B
       do ele_C=1,ele_count(libsupermesh_intersect_elements_result)
         area_intersect(index)  = area_intersect(index) + abs(simplex_volume(libsupermesh_intersect_elements_result, ele_C))
       end do
 
       call deallocate(libsupermesh_intersect_elements_result)
    end do
  end do
  
  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)
       allocate(posA(positionsA%dim, ele_loc(positionsA, ele_A)))
       allocate(posB(positionsB%dim, ele_loc(positionsB, ele_B)))
       posA = ele_val(positionsA, ele_A)
       posB = ele_val(positionsB, ele_B)
       elementShape = ele_shape(positionsB, 1)
       dimB = positionsB%dim
       locB = ele_loc(positionsB, ele_B)

       t1 = MPI_Wtime();
       call libsupermesh_cintersector_set_input(posB, posA, dimB, locB)
       call libsupermesh_cintersector_drive
       call libsupermesh_cintersector_query(nonods, totele)
       call allocate(intersection_meshLibWM, nonods, totele, elementShape, "IntersectionMeshLibWM")
       intersection_meshLibWM%continuity = -1
       call allocate(libwm_result, dimB, intersection_meshLibWM, "IntersectionCoordinatesLibWM")
       if (nonods > 0) then
         call libsupermesh_cintersector_get_output(nonods, totele, dimB, loc, nodes_tmp, intersection_meshLibWM%ndglno)

         do i = 1, dimB
           libwm_result%val(i,:) = nodes_tmp((i - 1) * nonods + 1:i * nonods)
         end do
       end if
       t2 = MPI_Wtime();
       dtLibWM = dtLibWM + ( t2 -t1 )

       deallocate(posA)
       deallocate(posB)
       call deallocate(intersection_meshLibWM)

       index = (ele_A-1)*BUF_SIZE_B + ele_B
       do ele_C=1,ele_count(libwm_result)
         area_libwm(index) = area_libwm(index) + abs(simplex_volume(libwm_result, ele_C))
       end do
 
       call deallocate(libwm_result)
    end do
  end do
  
  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)
      triA%v = ele_val(positionsA, ele_A)
      triB%v = ele_val(positionsB, ele_B)
      t1 = MPI_Wtime();
      call libsupermesh_intersect_tris_libwm(triA%v, triB%v, nodesC, ndglnoC, n_trisC)
      t2 = MPI_Wtime();
      dtLibSuperMeshIntersectLibWM = dtLibSuperMeshIntersectLibWM + ( t2 -t1 )

      index = (ele_A-1)*BUF_SIZE_B + ele_B
      do ele_C=1,n_trisC
        area_intersect_libwm(index) = area_intersect_libwm(index) + triangle_area(nodesC(:, ndglnoC(:, ele_C)))
      end do
    
    end do
  end do

  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)
      
      triA%v = ele_val(positionsA, ele_A)
      triB%v = ele_val(positionsB, ele_B)
      t1 = MPI_Wtime();
      call libsupermesh_intersect_tris_dt(triA, triB, trisC, n_trisC)
      t2 = MPI_Wtime();
      dtLibSuperMeshTriIntersector = dtLibSuperMeshTriIntersector + ( t2 -t1 )

      index = (ele_A-1)*BUF_SIZE_B + ele_B
      do ele_C=1,n_trisC
          area_fort(index) = area_fort(index) + triangle_area(trisC(ele_C)%v)
      end do
    end do
  end do
  
  write (*,*) "benchmark_tri_intersector: dtLibWM:",dtLibWM, &
      ", libSuperMeshIntersectLibWM:",dtLibSuperMeshIntersectLibWM, &
      ", dtLibSuperMesh::intersect_elements:",dtLibSuperMeshIntersectElements, &
      ", dtLibSuperMesh::tri:",dtLibSuperMeshTriIntersector,"."
  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)
      index = (ele_A-1)*BUF_SIZE_B + ele_B
      fail1 = (area_intersect(index) .fne. area_fort(index))
      fail2 = (area_libwm(index) .fne. area_intersect_libwm(index))
      
      if ( fail1 .OR. fail2 ) then
        write (*,*) "benchmark_tri_intersector: index:",index, &
           ", area_libwm:",area_libwm(index), &
           ", area_intersect_libwm:",area_intersect_libwm(index), &
           ", area_fort:",area_fort(index),&
           ", area_intersect:",area_intersect(index),"."
        write (*,*) "benchmark_tri_intersector: index:",index, &
           ", ele_A:",ele_A,", triA:",ele_val(positionsA, ele_A),&
           ", ele_B:",ele_B,", tri_B:",ele_val(positionsB, ele_B),"."
        totalFail = .TRUE.
        call report_test("[benchmark_tri_intersector areas]", fail1 .OR. fail2, .false., "Should give the same areas of intersection")
      end if
    end do
  end do
  
  call report_test("[benchmark_tri_intersector areas]", totalFail, .false., "Should give the same areas of intersection")
  
  call deallocate(positionsA)
  call deallocate(positionsB)
  
  call libsupermesh_cintersection_finder_reset(ntests)

contains

  pure function triangle_area(tri) result(area)
    real, dimension(2, 3), intent(in) :: tri

    real :: area
    real, dimension(2) :: u, v

    u = tri(:, 3) - tri(:, 1)
    v = tri(:, 2) - tri(:, 1)

    area = 0.5 * abs(u(2) * v(1) - u(1) * v(2))

  end function triangle_area

end subroutine benchmark_tri_intersector