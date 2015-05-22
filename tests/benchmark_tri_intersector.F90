#define BUF_SIZE_A 200
#define BUF_SIZE_B 200
subroutine benchmark_tri_intersector

  use libsupermesh_construction
  use libsupermesh_fields_allocates
  use libsupermesh_fields_base
  use libsupermesh_read_triangle
  use libsupermesh_tri_intersection_module
  use libsupermesh_unittest_tools
  
  implicit none

  interface mpi_wtime
    function mpi_wtime()
      implicit none
      real :: mpi_wtime
    end function mpi_wtime
  end interface mpi_wtime

  type(vector_field) :: positionsA, positionsB
  integer :: ele_A, ele_B, ele_C
  real, dimension(BUF_SIZE_A * BUF_SIZE_B) :: area_A_intersect_libwm,area_B_fort,area_C_fort_public,area_D_libwm,area_E_intersect_elements
  logical :: fail, totalFail = .FALSE.
  
  type(tri_type) :: triA, triB
  type(tri_type), dimension(8) :: trisC
  real, dimension(2, 3, tri_buf_size) ::  trisC_real
  integer :: ntests, n_trisC
  real, dimension(2, tri_buf_size) :: nodesC
  integer, dimension(3, tri_buf_size) :: ndglnoC
  
  type(vector_field) :: libsupermesh_intersect_elements_result, libwm_result
  real :: t1,t2,dt_A_area_intersect_libwm,dt_B_LibSuperMeshTriIntersector,dt_C_area_fort_public,dt_D_area_libwm,dt_E_area_intersect_elements
  real, dimension(:, :), allocatable :: posB, posA
  type(mesh_type) :: intersection_meshLibWM
  REAL :: num

  integer :: index, ele, dimB, i, locB, loc, elementCount, n_count, nodeCount
  type(element_type) :: elementShape
  
  real, dimension(:,:), allocatable :: positions_B_lib_val
  integer, dimension(3) :: counters
  
  integer :: nonods, totele
  
  ! I hope this is big enough ...
  real, dimension(1024) :: nodes_tmp

  call libsupermesh_intersector_set_dimension(2)
  call libsupermesh_intersector_set_exactness(.false.)
  
  open (unit = 20, file = "plcC_temp.node")
  open (unit = 21, file = "plcC_temp.ele")
  open (unit = 30, file = "plcD_temp.node")
  open (unit = 31, file = "plcD_temp.ele")

  i = 20
  CALL RANDOM_SEED(i)
  
  write (20,*) "",BUF_SIZE_A*3," 2 0 0"
  write (21,*) "",BUF_SIZE_A," 3 0"
  do ele=1,BUF_SIZE_A*3
    CALL RANDOM_NUMBER(num)
    num = num * 4.0
    if(num < 1.0) then
      write(20, *) " ", ele, " ", 1.0, num, " 0"
    else if(num < 2.0) then
      write(20, *) " ", ele, " ", 1.0 - (num - 1.0), 1.0, " 0"
    else if(num < 3.0) then
      write(20, *) " ", ele, " ", 0.0, 1.0 - (num - 2.0), " 0"
    else
      write(20, *) " ", ele, " ", num - 3.0, 0.0, " 0"
    end if
  end do
  do ele=1,BUF_SIZE_A
    write (21,*) " ",ele," ",(ele*3)-2," ", (ele*3)-1," ",(ele*3)
  end do
  
  write (30,*) "",BUF_SIZE_B*3," 2 0 0"
  write (31,*) "",BUF_SIZE_B," 3 0"
  do ele=1,BUF_SIZE_B*3
    CALL RANDOM_NUMBER(num)
    num = num * 4.0
    if(num < 1.0) then
      write(30, *) " ", ele, " ", 1.0, num, " 0"
    else if(num < 2.0) then
      write(30, *) " ", ele, " ", 1.0 - (num - 1.0), 1.0, " 0"
    else if(num < 3.0) then
      write(30, *) " ", ele, " ", 0.0, 1.0 - (num - 2.0), " 0"
    else
      write(30, *) " ", ele, " ", num - 3.0, 0.0, " 0"
    end if
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
  
  dt_A_area_intersect_libwm = 0.0
  dt_B_LibSuperMeshTriIntersector = 0.0
  dt_C_area_fort_public = 0.0
  dt_D_area_libwm = 0.0
  dt_E_area_intersect_elements = 0.0
  area_A_intersect_libwm = 0.0
  area_B_fort = 0.0
  area_C_fort_public = 0.0
  area_D_libwm = 0.0
  area_E_intersect_elements = 0.0

  counters = 0
  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)  
      triA%v = ele_val(positionsA, ele_A)
      triB%v = ele_val(positionsB, ele_B)

      ! A. Use libWM without creating temporary vector fields.
      t1 = MPI_Wtime();
      call libsupermesh_intersect_tris_libwm(triA%v, triB%v, nodesC, ndglnoC, n_trisC)
      t2 = MPI_Wtime();
      dt_A_area_intersect_libwm = dt_A_area_intersect_libwm + ( t2 -t1 )

      index = (ele_A-1)*BUF_SIZE_B + ele_B
      do ele_C=1,n_trisC
        area_A_intersect_libwm(index) = area_A_intersect_libwm(index) + triangle_area(nodesC(:, ndglnoC(:, ele_C)))
      end do


      ! B. Use the libSuperMesh internal triangle intersector (using derived types as input)
      t1 = MPI_Wtime();
      call libsupermesh_intersect_tris_dt(triA, triB, trisC, n_trisC)
      t2 = MPI_Wtime();
      dt_B_LibSuperMeshTriIntersector = dt_B_LibSuperMeshTriIntersector + ( t2 -t1 )

      do ele_C=1,n_trisC
          area_B_fort(index) = area_B_fort(index) + triangle_area(trisC(ele_C)%v)
      end do


      ! C. Use the libSuperMesh internal triangle intersector (using only reals as input)
      t1 = MPI_Wtime();
      call libsupermesh_intersect_tris_dt_public(triA%v, triB%v, trisC_real, n_trisC)
      t2 = MPI_Wtime();
      dt_C_area_fort_public = dt_C_area_fort_public + ( t2 -t1 )

      do ele_C=1,n_trisC
        area_C_fort_public(index) = area_C_fort_public(index) + triangle_area(trisC_real(:,:,ele_C))
      end do

      if(n_trisC == 0) then
        counters(1) = counters(1) + 1
      else if(n_trisC == 1) then
        counters(3) = counters(3) + 1
      else
        counters(2) = counters(2) + 1
      end if

      ! D. Use libWM directly and create temporary vector field
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
       dt_D_area_libwm = dt_D_area_libwm + ( t2 -t1 )

       deallocate(posA)
       deallocate(posB)
       call deallocate(intersection_meshLibWM)

       do ele_C=1,ele_count(libwm_result)
         area_D_libwm(index) = area_D_libwm(index) + abs(simplex_volume(libwm_result, ele_C))
       end do
       call deallocate(libwm_result)


      ! E. Use libWM directly and create temporary vector field      
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
       dt_E_area_intersect_elements = dt_E_area_intersect_elements + ( t2 -t1 )

       deallocate(positions_B_lib_val)
       deallocate(posB)

       do ele_C=1,ele_count(libsupermesh_intersect_elements_result)
         area_E_intersect_elements(index)  = area_E_intersect_elements(index) + abs(simplex_volume(libsupermesh_intersect_elements_result, ele_C))
       end do
 
       call deallocate(libsupermesh_intersect_elements_result)
    end do
  end do  

  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)
      index = (ele_A-1)*BUF_SIZE_B + ele_B
      fail = (area_A_intersect_libwm(index) .fne. area_B_fort(index)) &
         .OR. (area_E_intersect_elements(index) .fne. area_C_fort_public(index) ) &
         .OR. (area_E_intersect_elements(index) .fne. area_B_fort(index)) &
         .OR. (area_A_intersect_libwm(index) .fne. area_D_libwm(index)) &
         .OR. (area_A_intersect_libwm(index) .fne. area_C_fort_public(index) )
     
      if ( fail ) then
        write (*,*) "benchmark_tri_intersector: index:",index, &
           & ", area_A_intersect_libwm:",area_A_intersect_libwm(index), &
           & ", area_B_fort:",area_B_fort(index), &
           & ", area_C_fort_public:",area_C_fort_public(index), &
           & ", area_D_libwm:",area_D_libwm(index), &
           & ", area_E_intersect_elements:",area_E_intersect_elements(index),"."
        write (*,*) "benchmark_tri_intersector: index:",index, &
           & ", ele_A:",ele_A,", triA:",ele_val(positionsA, ele_A),&
           & ", ele_B:",ele_B,", tri_B:",ele_val(positionsB, ele_B),"."
        totalFail = .TRUE.
        call report_test("[benchmark_tri_intersector areas]", fail, .false., "Should give the same areas of intersection")
      end if
    end do
  end do

  write(*, *) "counters:",counters
  write(*, *) "dt_A_area_intersect_libwm      :",dt_A_area_intersect_libwm
  write(*, *) "dt_B_LibSuperMeshTriIntersector:",dt_B_LibSuperMeshTriIntersector
  write(*, *) "dt_C_area_fort_public          :",dt_C_area_fort_public
  write(*, *) "dt_D_area_libwm                :",dt_D_area_libwm
  write(*, *) "dt_E_area_intersect_elements   :",dt_E_area_intersect_elements
      
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
