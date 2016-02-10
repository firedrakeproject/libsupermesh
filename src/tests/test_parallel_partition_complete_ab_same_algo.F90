subroutine test_parallel_partition_complete_ab_same_algo

  use iso_c_binding, only : c_ptr, c_int8_t
  use iso_fortran_env, only : output_unit

  use libsupermesh_unittest_tools
  use libsupermesh_construction
  use libsupermesh_fields
  use libsupermesh_read_triangle
  use libsupermesh_tri_intersection
  use libsupermesh_unittest_tools
  use libsupermesh_intersection_finder
  use libsupermesh_read_halos
  use libsupermesh_halo_ownership
  use libsupermesh_parallel_supermesh, only : parallel_supermesh

  implicit none

#include <finclude/petsc.h90>

  character(len = 1024) :: buffer, hostname
  character(len = int(log10(real(huge(0)))) + 1) :: rank_character, nprocs_character
  integer :: ele_A, ele_B, ele_C, i, ierr, nprocs, n_trisC, rank, serial_ele_A, serial_ele_B, dp_extent, int_extent, test_parallel_ele_B, nintersections
  integer, parameter :: dim = 2, root = 0
  integer, dimension(:), allocatable :: ele_ownerA, ele_ownerB
  type(halo_type) :: halo
  type(intersections), dimension(:), allocatable :: map_AB
  type(tri_type) :: tri_A, tri_B
  type(tri_type), dimension(tri_buf_size) :: trisC
  type(vector_field) :: positionsA, positionsB
  real, parameter :: tol = 1.0e3 * epsilon(0.0)
  real :: t0, serial_time, parallel_time, serial_read_time, parallel_read_time
  logical :: fail

  real :: area_parallel, area_serial, integral_parallel, integral_serial
  real, dimension(:), allocatable, target :: valsB

  ! Data
  integer(kind = c_int8_t), dimension(:), allocatable :: data
  integer, dimension(2) :: data_b_header
  real, dimension(:), allocatable, target :: data_b_data

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
    positionsA = read_triangle_files("data/square_0_5"//"_"//trim(nprocs_character), dim)
    positionsB = read_triangle_files("data/square_0_9"//"_"//trim(nprocs_character), dim)
    serial_ele_A = ele_count(positionsA)
    serial_ele_B = ele_count(positionsB)
    serial_read_time = mpi_wtime() - t0

    t0 = mpi_wtime()
    call rtree_intersection_finder_set_input(positionsA%val, reshape(positionsA%mesh%ndglno, (/ele_loc(positionsA, 1), ele_count(positionsA)/)))

    allocate(valsB(serial_ele_B))
    do ele_B = 1, serial_ele_B
      valsB(ele_B) = sum(positionsB%val(1, ele_nodes(positionsB, ele_B))) / 3.0
    end do
    area_serial = 0.0
    integral_serial = 0.0

    do ele_B = 1, ele_count(positionsB)
      tri_B%v = ele_val(positionsB, ele_B)
      call rtree_intersection_finder_find(tri_B%v)
      call rtree_intersection_finder_query_output(nintersections)
      do i = 1, nintersections
        call rtree_intersection_finder_get_output(ele_A, i)
        tri_A%v = ele_val(positionsA, ele_A)

        call intersect_tris(tri_A, tri_B, trisC, n_trisC)

        do ele_C = 1, n_trisC
          area_serial = area_serial + triangle_area(trisC(ele_C)%v)
          integral_serial = integral_serial + valsB(ele_B) * triangle_area(trisC(ele_C)%v)
        end do
      end do
    end do

    deallocate(valsB)
    call deallocate(positionsA)
    call deallocate(positionsB)

    serial_time = mpi_wtime() - t0
  end if

  call MPI_Barrier(MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)

  t0 = mpi_wtime()
  positionsA = read_triangle_files(trim("data/square_0_5_")//trim(nprocs_character)//"_"//trim(rank_character), dim)
  call read_halo("data/square_0_5"//"_"//trim(nprocs_character), halo, level = 2)
  allocate(ele_ownerA(ele_count(positionsA)))
  call element_ownership(node_count(positionsA), reshape(positionsA%mesh%ndglno, (/ele_loc(positionsA, 1), ele_count(positionsA)/)), halo, ele_ownerA)
  call deallocate(halo)

  positionsB = read_triangle_files(trim("data/square_0_9_")//trim(nprocs_character)//"_"//trim(rank_character), dim)
  call read_halo("data/square_0_9"//"_"//trim(nprocs_character), halo, level = 2)
  allocate(ele_ownerB(ele_count(positionsB)))
  call element_ownership(node_count(positionsB), reshape(positionsB%mesh%ndglno, (/ele_loc(positionsB, 1), ele_count(positionsB)/)), halo, ele_ownerB)
  test_parallel_ele_B = count(ele_ownerB == rank)
  call deallocate(halo)
  parallel_read_time = mpi_wtime() - t0

  t0 = mpi_wtime()
  allocate(valsB(test_parallel_ele_B))
  do ele_B = 1, test_parallel_ele_B
    valsB(ele_B) = sum(positionsB%val(1, ele_nodes(positionsB, ele_B))) / 3.0
  end do
  area_parallel = 0.0
  integral_parallel = 0.0
  call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, dp_extent, ierr); CHKERRQ(ierr)
  call MPI_TYPE_EXTENT(MPI_INTEGER, int_extent, ierr); CHKERRQ(ierr)

  call parallel_supermesh(positionsA%val, & 
           &  reshape(positionsA%mesh%ndglno, (/ele_loc(positionsA, 1), ele_count(positionsA)/)), &
           &  ele_ownerA,          &
           &  positionsB%val,            &
           &  reshape(positionsB%mesh%ndglno, (/ele_loc(positionsB, 1), ele_count(positionsB)/)), &
           &  ele_ownerB,          &
           &  local_donor_ele_data, local_unpack_data_b, local_intersection_calculation)
  parallel_time = mpi_wtime() - t0

  call MPI_Allreduce(MPI_IN_PLACE, area_parallel, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  call MPI_Allreduce(MPI_IN_PLACE, integral_parallel, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  call MPI_Allreduce(MPI_IN_PLACE, parallel_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  call MPI_Allreduce(MPI_IN_PLACE, parallel_read_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)

  call deallocate(positionsA)
  call deallocate(positionsB)

  call MPI_Barrier(MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)

  if(rank == root) then
    write(output_unit, "(a,f19.15)") "Time, serial         =", serial_time
    write(output_unit, "(a,f19.15)") "(MAX) Time, parallel =", parallel_time
    write(output_unit, "(a)") ""
    write(output_unit, "(a,f19.15)") "Read Time, serial         =", serial_read_time
    write(output_unit, "(a,f19.15)") "(MAX) Read Time, parallel =", parallel_read_time
    write(output_unit, "(a)") ""
    write(output_unit, "(a,f19.15)") "Area, serial   =", area_serial
    write(output_unit, "(a,f19.15)") "Area, parallel =", area_parallel
    write(output_unit, "(a)") ""
    write(output_unit, "(a,f19.15)") "Integral, serial   =", integral_serial
    write(output_unit, "(a,f19.15)") "Integral, parallel =", integral_parallel

    fail = fnequals(area_parallel, area_serial, tol = tol)
    call report_test("[test_parallel_partition_complete_ab areas]", fail, .FALSE., "Should give the same areas of intersection")
    if (fail) then
      print "(a,e25.17e3,a,e25.17e3,a)", ": Area, serial   =", area_serial, &
                                       & ": Area, parallel =", area_parallel
    end if

    fail = fnequals(integral_parallel, integral_serial, tol = tol)
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
    real, dimension(:), pointer :: data_b
    
    if(local) then
      data_b => valsB
    else
      data_b => data_b_data
    end if

    nelements_c = size(positions_c, 3)
    do ele_c = 1, nelements_c
      area_parallel = area_parallel + triangle_area(positions_c(:, :, ele_c))
      integral_parallel = integral_parallel + data_b(ele_b) * triangle_area(positions_c(:, :, ele_c))
    end do

  end subroutine local_intersection_calculation

  subroutine write_parallel(msg)
    character(len = *), intent(in) :: msg

    integer :: i, ierr

    call MPI_Barrier(MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    flush(output_unit)
    do i = 1, nprocs
      if(i == rank + 1) then
        write(output_unit, "(a,a,a)", advance = "no") "Rank ", trim(rank_character), ": "
        write(output_unit, "(a)") trim(msg)
        flush(output_unit)
      end if
      call MPI_Barrier(MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
      flush(output_unit)
    end do

    flush(output_unit)

  end subroutine write_parallel

end subroutine test_parallel_partition_complete_ab_same_algo
