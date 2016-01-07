subroutine benchmark_serial_same_algo_3D

  use iso_c_binding, only : c_ptr, c_int8_t
  use iso_fortran_env, only : output_unit

  use libsupermesh_unittest_tools
  use libsupermesh_construction
  use libsupermesh_fields
  use libsupermesh_read_triangle
  use libsupermesh_tet_intersection_module
  use libsupermesh_tri_intersection_module, only : cintersection_finder_reset
  use libsupermesh_unittest_tools
  use libsupermesh_intersection_finder
  use libsupermesh_read_halos
  use libsupermesh_halo_ownership
  use libsupermesh_parallel_supermesh, only : parallel_supermesh

  implicit none

#include <finclude/petsc.h90>

  character(len = 1024) :: buffer, hostname
  character(len = int(log10(real(huge(0)))) + 1) :: rank_character, nprocs_character
  integer :: ele_A, ele_B, ele_C, i, ierr, nprocs, n_tetsC, rank, serial_ele_A, serial_ele_B, nintersections
  integer, parameter :: dim = 3, root = 0
  type(tet_type) :: tet_A, tet_B
  type(tet_type), dimension(tet_buf_size) :: tetsC
  type(vector_field) :: positionsA, positionsB
  real, parameter :: tol = 1.0e3 * epsilon(0.0)
  real :: t0, serial_time, serial_read_time
  logical :: fail

  real :: vols_parallel, vols_serial, integral_parallel, integral_serial
  real, dimension(:), allocatable :: valsB

  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr); CHKERRQ(ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr); CHKERRQ(ierr)

  write(rank_character, "(i0)") rank
  rank_character = adjustl(rank_character)
  write(nprocs_character, "(i0)") nprocs
  nprocs_character = adjustl(nprocs_character)
#ifdef __GNUC__
  call hostnm(hostname)
#else
  hostname = "unknown"
#endif
  write(buffer, "(a,a,a,i0,a)") "Running on '", trim(hostname), "' with ", nprocs, " processes"
  if ((rank == 0) .or. (mod(rank,10) == 0)) print *, trim(buffer)

  ! Serial test
  if (rank == 0) then
    t0 = mpi_wtime()
    positionsA = read_triangle_files("data/pyramid_0_5"//"_"//trim(nprocs_character), dim)
    positionsB = read_triangle_files("data/pyramid_0_9"//"_"//trim(nprocs_character), dim)
    serial_ele_A = ele_count(positionsA)
    serial_ele_B = ele_count(positionsB)
    serial_read_time = mpi_wtime() - t0

    t0 = mpi_wtime()
    call rtree_intersection_finder_set_input(positionsA%val, reshape(positionsA%mesh%ndglno, (/ele_loc(positionsA, 1), ele_count(positionsA)/)))

    allocate(valsB(serial_ele_B))
    do ele_B = 1, serial_ele_B
      valsB(ele_B) = sum(positionsB%val(1, ele_nodes(positionsB, ele_B))) / 4.0
    end do
    vols_serial = 0.0
    integral_serial = 0.0

    do ele_B = 1, ele_count(positionsB)
      tet_B%v = ele_val(positionsB, ele_B)
      call rtree_intersection_finder_find(tet_B%v)
      call rtree_intersection_finder_query_output(nintersections)
      do i = 1, nintersections
        call rtree_intersection_finder_get_output(ele_A, i)
        tet_A%v = ele_val(positionsA, ele_A)

        call intersect_tets(tet_A, tet_B, tetsC, n_tetsC)

        do ele_C = 1, n_tetsC
          vols_serial = vols_serial + tetrahedron_volume(tetsC(ele_C)%v)
          integral_serial = integral_serial + valsB(ele_B) * tetrahedron_volume(tetsC(ele_C)%v)
        end do
      end do
    end do

    deallocate(valsB)
    call deallocate(positionsA)
    call deallocate(positionsB)
    call cintersection_finder_reset(nintersections)

    serial_time = mpi_wtime() - t0
  end if

  call MPI_Allreduce(MPI_IN_PLACE, vols_parallel, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  call MPI_Allreduce(MPI_IN_PLACE, integral_parallel, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  vols_parallel = 1125.0
  integral_parallel = 5625.0

  if(rank == root) then
    write(output_unit, "(a,f19.15)") "Time, serial         =", serial_time
    write(output_unit, "(a)") ""
    write(output_unit, "(a,f19.15)") "Read Time, serial    =", serial_read_time
    write(output_unit, "(a)") ""
    write(output_unit, "(a,f19.14)") "Volume, serial       =", vols_serial
    write(output_unit, "(a,f19.14)") "Integral, serial     =", integral_serial

    fail = fnequals(vols_parallel, vols_serial, tol = tol)
    call report_test("[test_parallel_partition_complete_ab areas]", fail, .FALSE., "Should give the same areas of intersection")
    if (fail) then
      print "(a,e25.17e3,a,e25.17e3,a)", ": Volume, serial   =", vols_serial, &
                                       & ": Volume, parallel =", vols_parallel
    end if

    fail = fnequals(integral_parallel, integral_serial, tol = tol)
    call report_test("[test_parallel_partition_complete_ab integrals]", fail, .FALSE., "Should give the same values of integration")
    if (fail) then
      print "(a,e25.17e3,a,e25.17e3,a)", ": Integral, serial   =", integral_serial, &
                                       & ": Integral, parallel =", integral_parallel
    end if
  end if

end subroutine benchmark_serial_same_algo_3D