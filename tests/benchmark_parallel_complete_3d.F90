#include "libsupermesh_debug.h"

subroutine benchmark_parallel_complete_3D() bind(c)

  use iso_c_binding, only : c_ptr, c_int8_t
  use iso_fortran_env, only : output_unit

  use libsupermesh_debug
  use libsupermesh_unittest_tools
  use libsupermesh_supermesh
  use libsupermesh_read_triangle
  use libsupermesh_tri_intersection
  use libsupermesh_unittest_tools
  use libsupermesh_intersection_finder
  use libsupermesh_read_halos
  use libsupermesh_halo_ownership
  use libsupermesh_parallel_supermesh, only : parallel_supermesh, print_times

  implicit none

#include <mpif.h>

  character(len = 1024) :: buffer, hostname
  character(len = int(log10(real(huge(0)))) + 1) :: rank_character, nprocs_character
  integer :: ele_B, ierr, nprocs, rank, test_parallel_ele_B, dp_extent, int_extent
  integer, parameter :: dim = 3, root = 0
  integer, dimension(:), allocatable :: ele_ownerA, ele_ownerB
  type(halo_type) :: halo
  integer, dimension(:, :), allocatable :: enlist_a, enlist_b
  real, dimension(:, :), allocatable :: positions_a, positions_b

  real :: t0, serial_time, parallel_time, serial_read_time, parallel_read_time
  real :: parallel_time_tot_min, parallel_time_tot_max, parallel_time_tot_sum 
  logical :: fail

  real :: vols_parallel, vols_serial, integral_parallel, integral_serial
  real, dimension(:), allocatable :: valsB
#ifdef PROFILE
  real :: all_to_all_max, all_to_all_min, all_to_all_sum
  real :: point_to_point_max, point_to_point_min, point_to_point_sum
#endif

  ! Data
  integer(kind = c_int8_t), dimension(:), allocatable :: data
  integer, dimension(2) :: data_b_header
  real, dimension(:), allocatable :: data_b_data

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
  if ((rank == 0) .or. (mod(rank,20) == 0)) print *, trim(buffer)

  t0 = mpi_wtime()
  call read_node("data/pyramid_0_5_" // trim(nprocs_character) // "_" // trim(rank_character) // ".node", dim, positions_a)
  call read_ele("data/pyramid_0_5_" // trim(nprocs_character) // "_" // trim(rank_character) // ".ele", dim, enlist_a)
  call read_halo("data/pyramid_0_5_"//trim(nprocs_character), halo, level = 2)
  allocate(ele_ownerA(size(enlist_a, 2)))
  call element_ownership(size(positions_a, 2), enlist_a, halo, ele_ownerA)
  call deallocate(halo)

  call read_node("data/cube_0_5_" // trim(nprocs_character) // "_" // trim(rank_character) // ".node", dim, positions_b)
  call read_ele("data/cube_0_5_" // trim(nprocs_character) // "_" // trim(rank_character) // ".ele", dim, enlist_b)
  call read_halo("data/cube_0_5_"//trim(nprocs_character), halo, level = 2)
  allocate(ele_ownerB(size(enlist_b, 2)))
  call element_ownership(size(positions_b, 2), enlist_b, halo, ele_ownerB)
  test_parallel_ele_B = count(ele_ownerB == rank)
  call deallocate(halo)
  parallel_read_time = mpi_wtime() - t0

  if (rank == 0) write(output_unit, *) " A:", size(enlist_a, 2), ", B:", size(enlist_b, 2)

  t0 = mpi_wtime()
  allocate(valsB(test_parallel_ele_B))
  do ele_B = 1, test_parallel_ele_B
    valsB(ele_B) = sum(positions_b(1, enlist_b(:, ele_B))) / 4.0D0
  end do
  call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, dp_extent, ierr);  assert(ierr == MPI_SUCCESS)
  call MPI_TYPE_EXTENT(MPI_INTEGER, int_extent, ierr);  assert(ierr == MPI_SUCCESS)
  vols_parallel = 0.0D0
  integral_parallel = 0.0D0

  call parallel_supermesh(positions_a, & 
           &  enlist_a, &
           &  ele_ownerA,          &
           &  positions_b,            &
           &  enlist_b, &
           &  ele_ownerB,          &
           &  local_donor_ele_data, local_unpack_data_b, local_intersection_calculation)
  parallel_time = mpi_wtime() - t0

  call MPI_Allreduce(MPI_IN_PLACE, vols_parallel, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
  call MPI_Allreduce(MPI_IN_PLACE, integral_parallel, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
  call MPI_Allreduce(parallel_time, parallel_time_tot_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
  call MPI_Allreduce(parallel_time, parallel_time_tot_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
  call MPI_Allreduce(parallel_time, parallel_time_tot_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
  call MPI_Allreduce(MPI_IN_PLACE, parallel_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
  call MPI_Allreduce(MPI_IN_PLACE, parallel_read_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)

  deallocate(positions_a, enlist_a)
  deallocate(positions_b, enlist_b)

  vols_serial = 1000.0D0 / 3.0D0
  integral_serial = 0.166666074270734271D+004

  call print_times()

  if(rank == root) then
    write(output_unit, "(a)") ""
    write(output_unit, "(a,f20.10)") "Time, serial         = ", serial_time
    write(output_unit, "(a,f20.10)") "(MIN) Time, parallel = ", parallel_time_tot_min
    write(output_unit, "(a,f20.10)") "(MAX) Time, parallel = ", parallel_time_tot_max
    write(output_unit, "(a,f20.10)") "(SUM) Time, parallel = ", parallel_time_tot_sum
    write(output_unit, "(a,f20.10)") "(AVG) Time, parallel = ", parallel_time_tot_sum / nprocs
    write(output_unit, "(a)") ""
    write(output_unit, "(a,f19.10)") "Read Time, serial         = ", serial_read_time
    write(output_unit, "(a,f19.10)") "(MAX) Read Time, parallel = ", parallel_read_time
    write(output_unit, "(a)") ""
    write(output_unit, "(a,f19.14)") "Volume, parallel = ", vols_parallel
    write(output_unit, "(a)") ""
    write(output_unit, "(a,f19.14)") "Integral, parallel = ", integral_parallel
    write(output_unit, "(a)") ""

    fail = fnequals(vols_parallel, vols_serial)
    call report_test("[test_parallel_partition_complete_ab areas]", fail, .FALSE., "Should give the same areas of intersection")
    if (fail) then
      print "(a,e25.17e3,a,e25.17e3,a)", ": Volume, serial   =", vols_serial, &
                                       & ": Volume, parallel =", vols_parallel
    end if

    fail = fnequals(integral_parallel, integral_serial, tol = 1.0D3 * spacing(integral_serial))
    call report_test("[test_parallel_partition_complete_ab integrals]", fail, .FALSE., "Should give the same values of integration")
    if (fail) then
      print "(a,e25.17e3,a,e25.17e3,a)", ": Integral, serial   =", integral_serial, &
                                       & ": Integral, parallel =", integral_parallel
    end if
  end if

contains

  subroutine local_donor_ele_data(nodes, eles, data)
    integer, dimension(:), intent(in)                                :: nodes
    integer, dimension(:), intent(in)                                :: eles
    integer(kind = c_int8_t), dimension(:), allocatable, intent(out) :: data

    integer :: ierr, ndata, position
    real, dimension(:), allocatable :: ldata

    allocate(ldata(size(eles)))
    ldata = valsB(eles)

    ndata = 2 * int_extent + size(eles) * dp_extent
    allocate(data(ndata))
    position = 0
    call MPI_Pack((/size(eles), size(eles)/), 2, MPI_INTEGER, data, ndata, position, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
    call MPI_Pack(ldata, size(eles), MPI_DOUBLE_PRECISION, data, ndata, position, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)

    deallocate(ldata)

  end subroutine local_donor_ele_data

  subroutine local_unpack_data_b(nnodes_b, nelements_b, data_b)
    integer, intent(in) :: nnodes_b
    integer, intent(in) :: nelements_b
    integer(kind = c_int8_t), dimension(:), intent(in) :: data_b

    integer :: ierr, position

    call cleanup_data_b(data_b)

    position = 0
    call MPI_Unpack(data_b, size(data_b), position, data_b_header, 2, MPI_INTEGER, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
    allocate(data_b_data(data_b_header(2)))
    call MPI_Unpack(data_b, size(data_b), position, data_b_data, size(data_b_data), MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)

  end subroutine local_unpack_data_b

  subroutine cleanup_data_b(data_b)
    integer(kind = c_int8_t), dimension(:), intent(in) :: data_b

    if(allocated(data_b_data)) deallocate(data_b_data)

  end subroutine cleanup_data_b

  subroutine local_intersection_calculation(positions_a, positions_b, positions_c, nodes_b, ele_a, ele_b, local)
    ! dim x loc_a
    real, dimension(:, :), intent(in) :: positions_a
    ! dim x loc_b
    real, dimension(:, :), intent(in) :: positions_b
    ! dim x loc_c x nelements_c
    real, dimension(:, :, :), intent(in) :: positions_c
    ! loc_b
    integer, dimension(:), intent(in) :: nodes_b
    integer, intent(in)     :: ele_a
    integer, intent(in)     :: ele_b
    logical, intent(in)     :: local

    integer :: ele_c, nelements_c

    nelements_c = size(positions_c, 3)

    if (local) then
      do ele_c = 1, nelements_c
        vols_parallel = vols_parallel + tetrahedron_volume(positions_c(:, :, ele_c))
        integral_parallel = integral_parallel + valsB(ele_b) * tetrahedron_volume(positions_c(:, :, ele_c))
      end do
    else
      do ele_c = 1, nelements_c
        vols_parallel = vols_parallel + tetrahedron_volume(positions_c(:, :, ele_c))
        integral_parallel = integral_parallel + data_b_data(ele_b) * tetrahedron_volume(positions_c(:, :, ele_c))
      end do
    end if

  end subroutine local_intersection_calculation

end subroutine benchmark_parallel_complete_3D
