subroutine benchmark_parallel_complete_3D

  use iso_c_binding, only : c_ptr, c_int8_t
  use iso_fortran_env, only : output_unit

  use libsupermesh_unittest_tools
  use libsupermesh_construction
  use libsupermesh_fields
  use libsupermesh_read_triangle
  use libsupermesh_tri_intersection_module
  use libsupermesh_unittest_tools
  use libsupermesh_intersection_finder
  use libsupermesh_read_halos
  use libsupermesh_halo_ownership
  use libsupermesh_parallel_supermesh, only : parallel_supermesh

  implicit none

#include <finclude/petsc.h90>

  character(len = 1024) :: buffer, hostname
  character(len = int(log10(real(huge(0)))) + 1) :: rank_character, nprocs_character
  integer :: ele_B, ierr, nprocs, rank, test_parallel_ele_B, dp_extent, int_extent
  integer, parameter :: dim = 3, root = 0
  integer, dimension(:), allocatable :: ele_ownerA, ele_ownerB
  type(halo_type) :: halo
  type(vector_field) :: positionsA, positionsB
  real, parameter :: tol = 1.0e3 * epsilon(0.0)
  real :: t0, serial_time, parallel_time, serial_read_time, parallel_read_time
  logical :: fail

  real :: vols_parallel, vols_serial, integral_parallel, integral_serial
  real, dimension(:), allocatable :: valsB

  ! Data
  integer(kind = c_int8_t), dimension(:), allocatable :: data
  integer, dimension(2) :: data_b_header
  real, dimension(:), allocatable :: data_b_data

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

  t0 = mpi_wtime()
  positionsA = read_triangle_files(trim("data/pyramid_0_5_")//trim(nprocs_character)//"_"//trim(rank_character), dim)
  call read_halo("data/pyramid_0_5"//"_"//trim(nprocs_character), halo, level = 2)
  allocate(ele_ownerA(ele_count(positionsA)))
  call element_ownership(node_count(positionsA), reshape(positionsA%mesh%ndglno, (/ele_loc(positionsA, 1), ele_count(positionsA)/)), halo, ele_ownerA)
  call deallocate(halo)

  positionsB = read_triangle_files(trim("data/pyramid_0_9_")//trim(nprocs_character)//"_"//trim(rank_character), dim)
  call read_halo("data/pyramid_0_9"//"_"//trim(nprocs_character), halo, level = 2)
  allocate(ele_ownerB(ele_count(positionsB)))
  call element_ownership(node_count(positionsB), reshape(positionsB%mesh%ndglno, (/ele_loc(positionsB, 1), ele_count(positionsB)/)), halo, ele_ownerB)
  test_parallel_ele_B = count(ele_ownerB == rank)
  call deallocate(halo)
  parallel_read_time = mpi_wtime() - t0

  t0 = mpi_wtime()
  allocate(valsB(test_parallel_ele_B))
  do ele_B = 1, test_parallel_ele_B
    valsB(ele_B) = sum(positionsB%val(1, ele_nodes(positionsB, ele_B))) / 4.0
  end do
  call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, dp_extent, ierr);  CHKERRQ(ierr)
  call MPI_TYPE_EXTENT(MPI_INTEGER, int_extent, ierr);  CHKERRQ(ierr)
  vols_parallel = 0.0
  integral_parallel = 0.0

  call parallel_supermesh(positionsA%val, & 
           &  reshape(positionsA%mesh%ndglno, (/ele_loc(positionsA, 1), ele_count(positionsA)/)), &
           &  ele_ownerA,          &
           &  positionsB%val,            &
           &  reshape(positionsB%mesh%ndglno, (/ele_loc(positionsB, 1), ele_count(positionsB)/)), &
           &  ele_ownerB,          &
           &  local_donor_ele_data, local_unpack_data_b, local_intersection_calculation)
  parallel_time = mpi_wtime() - t0

  call MPI_Allreduce(MPI_IN_PLACE, vols_parallel, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  call MPI_Allreduce(MPI_IN_PLACE, integral_parallel, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  call MPI_Allreduce(MPI_IN_PLACE, parallel_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  call MPI_Allreduce(MPI_IN_PLACE, parallel_read_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)

  call deallocate(positionsA)
  call deallocate(positionsB)

  vols_serial = 1125.0
  integral_serial = 5625.0

  if(rank == root) then
    write(output_unit, "(a,f19.15)") "Time, serial         =", serial_time
    write(output_unit, "(a,f19.15)") "(MAX) Time, parallel =", parallel_time
    write(output_unit, "(a)") ""
    write(output_unit, "(a,f19.15)") "Read Time, serial         =", serial_read_time
    write(output_unit, "(a,f19.15)") "(MAX) Read Time, parallel =", parallel_read_time
    write(output_unit, "(a)") ""
    write(output_unit, "(a,f19.14)") "Volume, parallel =", vols_parallel
    write(output_unit, "(a)") ""
    write(output_unit, "(a,f19.14)") "Integral, parallel =", integral_parallel

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

contains

  subroutine local_donor_ele_data(eles, data)
    integer, dimension(:), intent(in)                   :: eles
    integer(kind = c_int8_t), dimension(:), allocatable :: data

    integer :: ierr, ndata, position
    real, dimension(:), allocatable :: ldata

    allocate(ldata(size(eles)))
    ldata = valsB(eles)

    ndata = 2 * int_extent + size(eles) * dp_extent
    allocate(data(ndata))
    position = 0
    call MPI_Pack((/size(eles), size(eles)/), 2, MPI_INTEGER, data, ndata, position, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    call MPI_Pack(ldata, size(eles), MPI_DOUBLE_PRECISION, data, ndata, position, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)

    deallocate(ldata)

  end subroutine local_donor_ele_data

  subroutine local_unpack_data_b(data_b)
    integer(kind = c_int8_t), dimension(:), intent(in) :: data_b

    integer :: ierr, position

    call cleanup_data_b(data_b)

    position = 0
    call MPI_Unpack(data_b, size(data_b), position, data_b_header, 2, MPI_INTEGER, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    allocate(data_b_data(data_b_header(2)))
    call MPI_Unpack(data_b, size(data_b), position, data_b_data, size(data_b_data), MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)

  end subroutine local_unpack_data_b

  subroutine cleanup_data_b(data_b)
    integer(kind = c_int8_t), dimension(:), intent(in) :: data_b

    if(allocated(data_b_data)) deallocate(data_b_data)

  end subroutine cleanup_data_b

  subroutine local_intersection_calculation(positions_c, ele_a, ele_b, local)
    ! dim x loc_c x nelements_c
    real, dimension(:, :, :), intent(in) :: positions_c
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
