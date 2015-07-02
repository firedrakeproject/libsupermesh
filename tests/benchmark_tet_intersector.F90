#define BUF_SIZE_A 1000
#define BUF_SIZE_B 1000
subroutine benchmark_tet_intersector

  use libsupermesh_construction
  use libsupermesh_fields_dummy
  use libsupermesh_read_triangle_2
  use libsupermesh_tet_intersection_module
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
  real, dimension(BUF_SIZE_A * BUF_SIZE_B) :: vol_A_libwm_intersect = 0.0, vol_B_fort = 0.0, &
  &  vol_C_fort_public = 0.0, vol_D_libwm = 0.0, vol_E_intersect_elements = 0.0, &
  &  vol_F_intersect_elements = 0.0, vol_G_intersect_elements = 0.0
  logical :: fail, totalFail = .FALSE.
  type(tet_type) :: tet_A, tet_B
  integer, dimension(8, BUF_SIZE_A * BUF_SIZE_B) :: tets_counters = 0
  integer, dimension(8) :: tets_totals = 0
  
  integer :: ele, i, index
  
  REAL, dimension(3) :: num
  double precision :: t1 = 0.0, t2 = 0.0, dt_A_vol_libwm_intersect = 0.0, dt_B_vol_fort = 0.0, &
  & dt_C_vol_fort_public = 0.0, dt_D_vol_libwm = 0.0, dt_E_vol_intersect_elements = 0.0, &
  & dt_F_vol_intersect_elements = 0.0, dt_G_vol_intersect_elements = 0.0
  integer :: nonods, totele
  
  real, dimension(3, 4, tet_buf_size) ::  tetsC_real
  integer :: n_tetsC
  real, dimension(3, tet_buf_size) :: nodesC
  integer, dimension(4, tet_buf_size) :: ndglnoC
  type(tet_type), dimension(tet_buf_size) :: tetsC
  type(mesh_type) :: intersection_mesh, new_mesh
  type(vector_field) :: intersection, intersect_elements_result
  integer, dimension(3) :: counters = 0
  
  ! I hope this is big enough ...
  real, dimension(1024) :: nodes_tmp

  integer, parameter :: dim = 3, loc = 4
  real, parameter :: tol = 1.0e3 * epsilon(0.0)

  call intersector_set_dimension(dim)
  call intersector_set_exactness(.false.)
  
  open (unit = 20, file = "plcC_temp.node")
  open (unit = 21, file = "plcC_temp.ele")
  open (unit = 30, file = "plcD_temp.node")
  open (unit = 31, file = "plcD_temp.ele")
  
  i = 20
  CALL RANDOM_SEED(i)
  
  write (20,*) "",BUF_SIZE_A*4," 3 0 0"
  write (21,*) "",BUF_SIZE_A," 4 0"
  do ele=1,BUF_SIZE_A*4
    CALL RANDOM_NUMBER(num)
    num = num * 10
    write (20,*) " ",ele," ",num
  end do
  do ele=1,BUF_SIZE_A
    write (21,*) " ",ele," ",(ele*4)-3," ",(ele*4)-2," ", (ele*4)-1," ",(ele*4)
  end do
  
  write (30,*) "",BUF_SIZE_B*4," 3 0 0"
  write (31,*) "",BUF_SIZE_B," 4 0"
  do ele=1,BUF_SIZE_B*4
    CALL RANDOM_NUMBER(num)
    num = num * 10
    write (30,*) " ",ele," ",num
  end do
  do ele=1,BUF_SIZE_B
    write (31,*) " ",ele," ",(ele*4)-3," ",(ele*4)-2," ", (ele*4)-1," ",(ele*4)
  end do
  
  close(20)
  close(21)
  close(30)
  close(31)
  
  positionsA = read_triangle_files("plcC_temp", dim)
  positionsB = read_triangle_files("plcD_temp", dim)
  
  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)
      ! A. Use libWM without creating temporary vector fields.
      tet_A%v = ele_val(positionsA, ele_A)
      tet_B%v = ele_val(positionsB, ele_B)

      t1 = MPI_Wtime();
      call intersect_tets_libwm(tet_A%v, tet_B%v, nodesC, ndglnoC, n_tetsC)
      t2 = MPI_Wtime();
      dt_A_vol_libwm_intersect = dt_A_vol_libwm_intersect + ( t2 -t1 )

      index = (ele_A-1)*BUF_SIZE_B + ele_B
      do ele_C=1,n_tetsC
        vol_A_libwm_intersect(index) = vol_A_libwm_intersect(index) + tetrahedron_volume(nodesC(:, ndglnoC(:, ele_C)))
      end do
      tets_counters(1, index) = n_tetsC
    end do
  end do


  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)
      ! B. Use the libSuperMesh internal triangle intersector (using derived types as input)
      tet_A%v = ele_val(positionsA, ele_A)
      tet_B%v = ele_val(positionsB, ele_B)
      
      t1 = MPI_Wtime();
      call intersect_tets_dt_public(tet_A, tet_B, tetsC, n_tetsC)
      t2 = MPI_Wtime();
      dt_B_vol_fort = dt_B_vol_fort + ( t2 -t1 )

      index = (ele_A-1)*BUF_SIZE_B + ele_B
      do ele_C=1,n_tetsC
        vol_B_fort(index) = vol_B_fort(index) + tetrahedron_volume(tetsC(ele_C)%v)
!        vol_A_libwm_intersect(index) = vol_B_fort(index)   ! HACK REMOVE
      end do
      tets_counters(2, index) = n_tetsC

      if(n_tetsC == 0) then
        counters(1) = counters(1) + 1
      else if(n_tetsC == 1) then
        counters(3) = counters(3) + 1
      else
        counters(2) = counters(2) + 1
      end if
    end do
  end do


  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)
      ! C. Use the libSuperMesh internal triangle intersector (using only reals as input)
      tet_A%v = ele_val(positionsA, ele_A)
      tet_B%v = ele_val(positionsB, ele_B)
      
      t1 = MPI_Wtime();
      call intersect_tets_dt_real(tet_A%v, tet_B%v, tetsC_real, n_tetsC)
      t2 = MPI_Wtime();
      dt_C_vol_fort_public = dt_C_vol_fort_public + ( t2 -t1 )

      index = (ele_A-1)*BUF_SIZE_B + ele_B
      do ele_C=1,n_tetsC
        vol_C_fort_public(index) = vol_C_fort_public(index) + tetrahedron_volume(tetsC_real(:,:,ele_C))
      end do
      tets_counters(3, index) = n_tetsC
    end do
  end do


  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)
      ! D. Use libWM directly and create temporary vector field
      tet_A%v = ele_val(positionsA, ele_A)
      tet_B%v = ele_val(positionsB, ele_B)

      t1 = MPI_Wtime();
      call libsupermesh_cintersector_set_input(tet_A%v, tet_B%v, dim, loc)
      call libsupermesh_cintersector_drive
      call libsupermesh_cintersector_query(nonods, totele)
      call allocate(intersection_mesh, dim, nonods, totele, loc)
      intersection_mesh%continuity = -1
      call allocate(intersection, dim, intersection_mesh)
      if (nonods > 0) then
        call libsupermesh_cintersector_get_output(nonods, totele, dim, loc, nodes_tmp, intersection_mesh%ndglno)
        do i = 1, dim
          intersection%val(i,:) = nodes_tmp((i - 1) * nonods + 1:i * nonods)
        end do
      end if
      t2 = MPI_Wtime();
      dt_D_vol_libwm = dt_D_vol_libwm + ( t2 -t1 )

      index = (ele_A-1)*BUF_SIZE_B + ele_B
      do ele_C=1,totele
        vol_D_libwm(index) = vol_D_libwm(index) + tetrahedron_volume(ele_val(intersection, ele_C))
      end do
      tets_counters(4, index) = ele_count(intersection)
      call deallocate(intersection_mesh)
      call deallocate(intersection)
    end do
  end do


  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)
      ! E. Use *original/old* intersect_elements function directly and create temporary vector field
      t1 = MPI_Wtime();
      intersect_elements_result = intersect_elements_old(ele_val(positionsA, ele_A), &
        ele_val(positionsB, ele_B), loc, dim)
      t2 = MPI_Wtime();
      dt_E_vol_intersect_elements = dt_E_vol_intersect_elements + ( t2 -t1 )

      index = (ele_A-1)*BUF_SIZE_B + ele_B
      do ele_C=1,ele_count(intersect_elements_result)
        vol_E_intersect_elements(index) = vol_E_intersect_elements(index) + tetrahedron_volume(ele_val(intersect_elements_result, ele_C))
      end do
      tets_counters(5, index) = ele_count(intersect_elements_result)
      call deallocate(intersect_elements_result)
    end do
  end do


  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)
      ! F. Use the new intersect_elements and do NOT create vector field
      tet_A%v = ele_val(positionsA, ele_A)
      tet_B%v = ele_val(positionsB, ele_B)

      t1 = MPI_Wtime();
      call intersect_elements(tet_A%v, tet_B%v, dim, n_tetsC, tetsC_real=tetsC_real)
      t2 = MPI_Wtime();
      dt_F_vol_intersect_elements = dt_F_vol_intersect_elements + ( t2 -t1 )

      index = (ele_A-1)*BUF_SIZE_B + ele_B
      do ele_C=1,n_tetsC
        vol_F_intersect_elements(index) = vol_F_intersect_elements(index) + tetrahedron_volume(tetsC_real(:,:,ele_C))
      end do
      tets_counters(6, index) = n_tetsC
    end do
  end do


  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)
      ! G. Use the new intersect_elements and DO create vector field
      tet_A%v = ele_val(positionsA, ele_A)
      tet_B%v = ele_val(positionsB, ele_B)

      t1 = MPI_Wtime();
      call intersect_elements(tet_A%v, tet_B%v, dim, n_tetsC, tetsC_real=tetsC_real)
      call allocate(new_mesh, dim, n_tetsC * loc, n_tetsC, loc)

      if ( n_tetsC > 0 ) then
        new_mesh%ndglno = (/ (i, i=1,loc * n_tetsC) /)
        new_mesh%continuity = -1
      end if

      call allocate(intersection, dim, new_mesh)
      if ( n_tetsC > 0 ) then
        do i = 1, n_tetsC
          call set(intersection, ele_nodes(intersection, i), tetsC_real(:,:,i))
        end do
      end if
      call deallocate(new_mesh)
      t2 = MPI_Wtime();
      dt_G_vol_intersect_elements = dt_G_vol_intersect_elements + ( t2 -t1 )

      index = (ele_A-1)*BUF_SIZE_B + ele_B
      do ele_C=1,ele_count(intersection)
        vol_G_intersect_elements(index) = vol_G_intersect_elements(index) + tetrahedron_volume(ele_val(intersection, ele_C))
      end do
      tets_counters(7, index) = ele_count(intersection)
      call deallocate(intersection)
    end do
  end do


  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)
      index = (ele_A-1)*BUF_SIZE_B + ele_B
      do i = 1, 8
        tets_totals(i) = tets_totals(i) + tets_counters(i, index)
      end do
      fail = (fnequals(vol_A_libwm_intersect(index), vol_B_fort(index), tol = tol) &
         .OR. fnequals(vol_E_intersect_elements(index), vol_C_fort_public(index), tol = tol) &
         .OR. fnequals(vol_E_intersect_elements(index), vol_B_fort(index), tol = tol) &
         .OR. fnequals(vol_A_libwm_intersect(index), vol_D_libwm(index), tol = tol) &
         .OR. fnequals(vol_A_libwm_intersect(index), vol_C_fort_public(index), tol = tol) &
         .OR. fnequals(vol_F_intersect_elements(index), vol_A_libwm_intersect(index), tol = tol) &
         .OR. fnequals(vol_F_intersect_elements(index), vol_D_libwm(index), tol = tol) &
         .OR. fnequals(vol_G_intersect_elements(index), vol_B_fort(index), tol = tol) &
         .OR. fnequals(vol_G_intersect_elements(index), vol_D_libwm(index), tol = tol))
      
      if ( fail ) then
        write (*,*) "benchmark_tet_intersector: index:",index, &
          & ", vol_A_libwm_intersect:",vol_A_libwm_intersect(index),&
          & ", vol_B_fort:",vol_B_fort(index),&
          & ", vol_C_fort_public:",vol_C_fort_public(index),", vol_D_libwm:",vol_D_libwm(index),&
          & ", vol_E_intersect_elements:",vol_E_intersect_elements(index), &
          & ", vol_F_intersect_elements:",vol_F_intersect_elements(index), &
          & ", vol_G_intersect_elements:",vol_G_intersect_elements(index),"."
        write (*,*) "benchmark_tet_intersector: index:",index,     &
          & ", ele_A:",ele_A,", tetA:",ele_val(positionsA, ele_A), &
          & ", ele_B:",ele_B,", tet_B:",ele_val(positionsB, ele_B),"."
        totalFail = .TRUE.
        call report_test("[benchmark_tet_intersector areas]", fail, .false., "Should give the same volumes of intersection")
      end if
    end do
  end do

  write(*, *) "counters:",counters
  write(*, *) "A. Intersect libwm                  :",dt_A_vol_libwm_intersect,", tets:",tets_totals(1)
  write(*, *) "B. LibSuperMeshTetIntersector       :",dt_B_vol_fort,", tets:",tets_totals(2)
  write(*, *) "C. LibSuperMeshTetIntersector (real):",dt_C_vol_fort_public,", tets:",tets_totals(3)
  write(*, *) "D. Calling directly libwm           :",dt_D_vol_libwm,", tets:",tets_totals(4)," (create vector field)."
  write(*, *) "E. Intersect Elements (original)    :",dt_E_vol_intersect_elements,", tets:",tets_totals(5)," (create vector field)."
  write(*, *) "F. Intersect Elements (new)         :",dt_F_vol_intersect_elements,", tets:",tets_totals(6)
  write(*, *) "G. Intersect Elements (new)         :",dt_G_vol_intersect_elements,", tets:",tets_totals(7)," (create vector field)."
  
  call report_test("[benchmark_tet_intersector areas]", totalFail, .false., "Should give the same volumes of intersection")
  
  call deallocate(positionsA)
  call deallocate(positionsB)
  
  call libsupermesh_cintersection_finder_reset(n_tetsC)
  call finalise_libsupermesh()

end subroutine benchmark_tet_intersector
