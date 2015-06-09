#define BUF_SIZE_A 1000
#define BUF_SIZE_B 1000
subroutine benchmark_tri_intersector

  use libsupermesh_construction
  use libsupermesh_fields
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
  real, dimension(BUF_SIZE_A * BUF_SIZE_B) :: area_A_intersect_libwm = 0.0, &
    &  area_B_fort = 0.0, area_C_fort_public = 0.0, area_D_libwm = 0.0, &
    &  area_E_intersect_elements = 0.0, area_F_intersect_elements = 0.0, &
    &  area_G_intersect_elements = 0.0, area_H_intersect_old = 0.0
  logical :: fail, totalFail = .FALSE.
  integer, dimension(8, BUF_SIZE_A * BUF_SIZE_B) :: triangle_counters
  integer, dimension(8) :: triangle_totals
  
  type(tri_type) :: triA, triB
  type(tri_type), dimension(tri_buf_size) :: trisC
  real, dimension(2, 3, tri_buf_size) ::  trisC_real
  integer :: ntests, n_trisC
  real, dimension(2, tri_buf_size) :: nodesC
  integer, dimension(3, tri_buf_size) :: ndglnoC
  
  type(vector_field) :: libsupermesh_intersect_elements_result, libwm_result, intersection
  real :: t1 = 0.0, t2 = 0.0, dt_A_area_intersect_libwm = 0.0, dt_B_LibSuperMeshTriIntersector = 0.0, &
    &  dt_C_area_fort_public = 0.0, dt_D_area_libwm = 0.0, dt_E_area_intersect_elements = 0.0, &
    &  dt_F_area_intersect_elements = 0.0, dt_G_area_intersect_elements = 0.0, &
    &  dt_H_area_intersect_old = 0.0
  real, dimension(:, :), allocatable :: posB, posA
  type(mesh_type) :: intersection_meshLibWM, intersection_meshIntersect_Elements, new_mesh
  REAL :: num

  integer :: index, ele, dimB, i, locB, loc, elementCount, n_count, nodeCount
  type(element_type) :: elementShape
  
  real, dimension(:,:), allocatable :: positions_B_lib_val
  integer, dimension(3) :: counters
  
  integer :: nonods, totele, n_trisC_pre
  
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
  
  counters = 0

  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)  
      ! A. Use libWM without creating temporary vector fields.
      triA%v = ele_val(positionsA, ele_A)
      triB%v = ele_val(positionsB, ele_B)

      t1 = MPI_Wtime();
      call libsupermesh_intersect_tris_libwm(triA%v, triB%v, nodesC, ndglnoC, n_trisC)
      t2 = MPI_Wtime();
      dt_A_area_intersect_libwm = dt_A_area_intersect_libwm + ( t2 -t1 )

      index = (ele_A-1)*BUF_SIZE_B + ele_B
      do ele_C=1,n_trisC
        area_A_intersect_libwm(index) = area_A_intersect_libwm(index) + triangle_area(nodesC(:, ndglnoC(:, ele_C)))
      end do
      triangle_counters(1, index) = n_trisC
    end do
  end do


  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB) 
      ! B. Use the libSuperMesh internal triangle intersector (using derived types as input)
      triA%v = ele_val(positionsA, ele_A)
      triB%v = ele_val(positionsB, ele_B)
      
      t1 = MPI_Wtime();
      call libsupermesh_intersect_tris_dt(triA, triB, trisC, n_trisC)
      t2 = MPI_Wtime();
      dt_B_LibSuperMeshTriIntersector = dt_B_LibSuperMeshTriIntersector + ( t2 -t1 )

      index = (ele_A-1)*BUF_SIZE_B + ele_B
      do ele_C=1,n_trisC
          area_B_fort(index) = area_B_fort(index) + triangle_area(trisC(ele_C)%v)
      end do
      triangle_counters(2, index) = n_trisC

       if(n_trisC == 0) then
         counters(1) = counters(1) + 1
       else if(n_trisC == 1) then
         counters(3) = counters(3) + 1
       else
         counters(2) = counters(2) + 1
       end if

      if (  triangle_counters(2, index) .ne.  triangle_counters(1, index) ) then
        if ( abs(area_A_intersect_libwm(index) - area_B_fort(index)) .le. epsilon(0.0) ) then
          cycle
        end if
        if ( abs(area_A_intersect_libwm(index) - area_B_fort(index)) .gt. epsilon(0.0) ) then
          write (*,*) "!!!!!!!! Different areas as well"
        end if
        write(*,*) "LibWM returned:", triangle_counters(1, index),", and tris returned:", triangle_counters(2, index),"."
        write(*,*) "Area: area_A_intersect_libwm:",area_A_intersect_libwm(index),", area_B_fort:",area_B_fort(index),"."
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
    end do
   end do


  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB) 
      ! C. Use the libSuperMesh internal triangle intersector (using only reals as input)
      triA%v = ele_val(positionsA, ele_A)
      triB%v = ele_val(positionsB, ele_B)

      t1 = MPI_Wtime();
      call libsupermesh_intersect_tris_dt_public(triA%v, triB%v, trisC_real, n_trisC)
      t2 = MPI_Wtime();
      dt_C_area_fort_public = dt_C_area_fort_public + ( t2 -t1 )

      index = (ele_A-1)*BUF_SIZE_B + ele_B
      do ele_C=1,n_trisC
        area_C_fort_public(index) = area_C_fort_public(index) + triangle_area(trisC_real(:,:,ele_C))
      end do
      triangle_counters(3, index) = n_trisC
    end do
  end do


  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB) 
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

      index = (ele_A-1)*BUF_SIZE_B + ele_B
      do ele_C=1,ele_count(libwm_result)
        area_D_libwm(index) = area_D_libwm(index) + abs(simplex_volume(libwm_result, ele_C))
      end do
      triangle_counters(4, index) = ele_count(libwm_result)
      call deallocate(libwm_result)
    end do
  end do


  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB) 
      ! E. Use *original/old* intersect_elements function directly and create temporary vector field      
      dimB = positionsB%dim
      n_count = 0
      select case(positionsB%field_type)
      case(FIELD_TYPE_NORMAL)
!        n_count = node_count(positionsB%mesh)
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
      libsupermesh_intersect_elements_result = libsupermesh_intersect_elements_old( &
        positions_B_lib_val, posB, loc, dimB, nodeCount, &
        elementShape%quadrature%vertices, elementShape%quadrature%dim, elementShape%quadrature%ngi, &
        elementShape%quadrature%degree, elementShape%loc, elementShape%dim, elementShape%degree)
      t2 = MPI_Wtime();
      dt_E_area_intersect_elements = dt_E_area_intersect_elements + ( t2 -t1 )

      index = (ele_A-1)*BUF_SIZE_B + ele_B
      do ele_C=1,ele_count(libsupermesh_intersect_elements_result)
        area_E_intersect_elements(index)  = area_E_intersect_elements(index) + abs(simplex_volume(libsupermesh_intersect_elements_result, ele_C))
      end do
      triangle_counters(5, index) = ele_count(libsupermesh_intersect_elements_result)
      call deallocate(libsupermesh_intersect_elements_result)
      deallocate(posB)
      deallocate(positions_B_lib_val)
    end do
  end do


  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB) 
      ! F. Use the new libsupermesh_intersect_elements and do NOT create vector field 
      triA%v = ele_val(positionsA, ele_A)
      triB%v = ele_val(positionsB, ele_B)

      t1 = MPI_Wtime();
      call libsupermesh_intersect_elements(triA%v, triB%v, positionsA%dim, n_trisC, trisC_real=trisC_real)
      t2 = MPI_Wtime();
      dt_F_area_intersect_elements = dt_F_area_intersect_elements + ( t2 -t1 )

      index = (ele_A-1)*BUF_SIZE_B + ele_B
      do ele_C=1,n_trisC
        area_F_intersect_elements(index) = area_F_intersect_elements(index) + triangle_area(trisC_real(:,:,ele_C))
      end do
      triangle_counters(6, index) = n_trisC
    end do
  end do


  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB) 
      ! G. Use the new libsupermesh_intersect_elements and DO create vector field 
      triA%v = ele_val(positionsA, ele_A)
      triB%v = ele_val(positionsB, ele_B)

      t1 = MPI_Wtime();
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
      t2 = MPI_Wtime();
      dt_G_area_intersect_elements = dt_G_area_intersect_elements + ( t2 -t1 )

      index = (ele_A-1)*BUF_SIZE_B + ele_B
      do ele_C=1,ele_count(intersection)
        area_G_intersect_elements(index) = area_G_intersect_elements(index) + abs(simplex_volume(intersection, ele_C))
      end do
      triangle_counters(7, index) = ele_count(intersection)
      call deallocate(intersection)
    end do
  end do


  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB) 
      ! H. Use the *OLD* libSuperMesh internal triangle intersector (using derived types as input)
      triA%v = ele_val(positionsA, ele_A)
      triB%v = ele_val(positionsB, ele_B)

      t1 = MPI_Wtime();
      call libsupermesh_intersect_tris_dt_old(triA, triB, trisC, n_trisC)
      t2 = MPI_Wtime();
      dt_H_area_intersect_old = dt_H_area_intersect_old + ( t2 -t1 )

      index = (ele_A-1)*BUF_SIZE_B + ele_B
      do ele_C=1,n_trisC
          area_H_intersect_old(index) = area_H_intersect_old(index) + triangle_area(trisC(ele_C)%v)
      end do
      triangle_counters(8, index) = n_trisC
    end do
   end do


  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)
      index = (ele_A-1)*BUF_SIZE_B + ele_B
      do i = 1, 8
        triangle_totals(i) = triangle_totals(i) + triangle_counters(i, index)
      end do
      fail = (area_A_intersect_libwm(index) .fne. area_B_fort(index)) &
         .OR. (area_E_intersect_elements(index) .fne. area_C_fort_public(index) ) &
         .OR. (area_E_intersect_elements(index) .fne. area_B_fort(index)) &
         .OR. (area_A_intersect_libwm(index) .fne. area_D_libwm(index)) &
         .OR. (area_A_intersect_libwm(index) .fne. area_C_fort_public(index) ) &
         .OR. (area_F_intersect_elements(index) .fne. area_A_intersect_libwm(index) ) &
         .OR. (area_F_intersect_elements(index) .fne. area_D_libwm(index) ) &
         .OR. (area_G_intersect_elements(index) .fne. area_B_fort(index)) &
         .OR. (area_G_intersect_elements(index) .fne. area_D_libwm(index)) & 
         .OR. (area_H_intersect_old(index) .fne. area_A_intersect_libwm(index)) &
         .OR. (area_H_intersect_old(index) .fne. area_C_fort_public(index)) 
     fail = .false.
      if ( fail ) then
        write (*,*) "benchmark_tri_intersector: index:",index, &
           & ", area_A_intersect_libwm:",area_A_intersect_libwm(index), &
           & ", area_B_fort:",area_B_fort(index), &
           & ", area_C_fort_public:",area_C_fort_public(index), &
           & ", area_D_libwm:",area_D_libwm(index), &
           & ", area_E_intersect_elements:",area_E_intersect_elements(index), &
           & ", area_F_intersect_elements:",area_F_intersect_elements(index), & 
           & ", area_G_intersect_elements:",area_G_intersect_elements(index), &
           & ", area_H_intersect_elements:",area_H_intersect_old(index),"."
        write (*,*) "benchmark_tri_intersector: index:",index, &
           & ", ele_A:",ele_A,", triA:",ele_val(positionsA, ele_A),&
           & ", ele_B:",ele_B,", tri_B:",ele_val(positionsB, ele_B),"."
        totalFail = .TRUE.
        call report_test("[benchmark_tri_intersector areas]", fail, .false., "Should give the same areas of intersection")
        call exit(1)
      end if
    end do
  end do

  write(*, *) "counters:",counters
  write(*, *) "A. Intersect libwm                   :",dt_A_area_intersect_libwm,", triangles:",triangle_totals(1)
  write(*, *) "B. LibSuperMeshTriIntersector        :",dt_B_LibSuperMeshTriIntersector,", triangles:",triangle_totals(2)
  write(*, *) "C. LibSuperMeshTriIntersector (real) :",dt_C_area_fort_public,", triangles:",triangle_totals(3)
  write(*, *) "D. Calling directly libwm            :",dt_D_area_libwm,", triangles:",triangle_totals(4)," (create vector field)."
  write(*, *) "E. Intersect Elements (original)     :",dt_E_area_intersect_elements,", triangles:",triangle_totals(5)," (create vector field)."
  write(*, *) "F. Intersect Elements (new)          :",dt_F_area_intersect_elements,", triangles:",triangle_totals(6)
  write(*, *) "G. Intersect Elements                :",dt_G_area_intersect_elements,", triangles:",triangle_totals(7)," (create vector field)."
  write(*, *) "H. LibSuperMeshTriIntersector (old)  :",dt_H_area_intersect_old,", triangles:",triangle_totals(8)
      
  call report_test("[benchmark_tri_intersector areas]", totalFail, .false., "Should give the same areas of intersection")
  
  call deallocate(positionsA)
  call deallocate(positionsB)
  
  call libsupermesh_cintersection_finder_reset(ntests)
  call finalise_libsupermesh()

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
