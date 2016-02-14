#include "libsupermesh_debug.h"

subroutine benchmark_serial_same_algo_3D() bind(c)

  use iso_c_binding, only : c_ptr, c_int8_t
  use iso_fortran_env, only : output_unit

  use libsupermesh_debug
  use libsupermesh_unittest_tools
  use libsupermesh_supermesh
  use libsupermesh_read_triangle
  use libsupermesh_tet_intersection
  use libsupermesh_unittest_tools
  use libsupermesh_intersection_finder
  use libsupermesh_read_halos
  use libsupermesh_halo_ownership

  implicit none

#include <mpif.h>

  character(len = 1024) :: buffer, hostname
  character(len = int(log10(real(huge(0)))) + 1) :: rank_character, nprocs_character
  integer :: ele_A, ele_B, ele_C, i, ierr, nprocs, n_tetsC, rank, serial_ele_A, serial_ele_B, nintersections
  integer, parameter :: dim = 3, root = 0
  type(tet_type) :: tet_A, tet_B
  type(tet_type), dimension(tet_buf_size) :: tetsC
  integer, dimension(:, :), allocatable :: enlist_a, enlist_b
  real, dimension(:, :), allocatable :: positions_a, positions_b
  real :: t0, serial_time, serial_read_time
  logical :: fail

  real :: vols_parallel, vols_serial, integral_parallel, integral_serial
  real, dimension(:), allocatable :: valsB

  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr);  assert(ierr == MPI_SUCCESS)
  CALL MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr);  assert(ierr == MPI_SUCCESS)

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
    call read_node("data/pyramid_0_5_" // trim(nprocs_character) // ".node", dim, positions_a)
    call read_ele("data/pyramid_0_5_" // trim(nprocs_character) // ".ele", dim, enlist_a)
    call read_node("data/pyramid_0_9_" // trim(nprocs_character) // ".node", dim, positions_b)
    call read_ele("data/pyramid_0_9_" // trim(nprocs_character) // ".ele", dim, enlist_b)
    serial_ele_A = size(enlist_a, 2)
    serial_ele_B = size(enlist_b, 2)
    serial_read_time = mpi_wtime() - t0

    if (rank == 0) write(output_unit, *) " A:", size(enlist_a, 2), ", B:", size(enlist_b, 2)

    t0 = mpi_wtime()
    call rtree_intersection_finder_set_input(positions_a, enlist_a)

    allocate(valsB(serial_ele_B))
    do ele_B = 1, serial_ele_B
      valsB(ele_B) = sum(positions_b(1, enlist_b(:, ele_B))) / 4.0D0
    end do
    vols_serial = 0.0D0
    integral_serial = 0.0D0

    do ele_B = 1, size(enlist_b, 2)
      tet_B%v = positions_b(:, enlist_b(:, ele_B))
      call rtree_intersection_finder_find(tet_B%v)
      call rtree_intersection_finder_query_output(nintersections)
      do i = 1, nintersections
        call rtree_intersection_finder_get_output(ele_A, i)
        tet_A%v = positions_a(:, enlist_a(:, ele_A))

        call intersect_tets(tet_A, tet_B, tetsC, n_tetsC)

        do ele_C = 1, n_tetsC
          vols_serial = vols_serial + tetrahedron_volume(tetsC(ele_C)%v)
          integral_serial = integral_serial + valsB(ele_B) * tetrahedron_volume(tetsC(ele_C)%v)
        end do
      end do
    end do

    deallocate(valsB)
    deallocate(positions_a, enlist_a)
    deallocate(positions_b, enlist_b)

    serial_time = mpi_wtime() - t0
  end if

  call MPI_Allreduce(MPI_IN_PLACE, vols_parallel, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
  call MPI_Allreduce(MPI_IN_PLACE, integral_parallel, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
  vols_parallel = 1125.0D0
  integral_parallel = 5625.0D0

  if(rank == root) then
    write(output_unit, "(a,f19.10)") "Time, serial         =", serial_time
    write(output_unit, "(a)") ""
    write(output_unit, "(a,f19.10)") "Read Time, serial    =", serial_read_time
    write(output_unit, "(a)") ""
    write(output_unit, "(a,f19.14)") "Volume, serial       =", vols_serial
    write(output_unit, "(a,f19.14)") "Integral, serial     =", integral_serial

    fail = fnequals(vols_parallel, vols_serial)
    call report_test("[test_parallel_partition_complete_ab areas]", fail, .FALSE., "Should give the same areas of intersection")
    if (fail) then
      print "(a,e25.17e3,a,e25.17e3,a)", ": Volume, serial   =", vols_serial, &
                                       & ": Volume, parallel =", vols_parallel
    end if

    fail = fnequals(integral_parallel, integral_serial)
    call report_test("[test_parallel_partition_complete_ab integrals]", fail, .FALSE., "Should give the same values of integration")
    if (fail) then
      print "(a,e25.17e3,a,e25.17e3,a)", ": Integral, serial   =", integral_serial, &
                                       & ": Integral, parallel =", integral_parallel
    end if
  end if

end subroutine benchmark_serial_same_algo_3D
