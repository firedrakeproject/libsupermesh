subroutine test_parallel_partition_complete_ab

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
  integer :: ele_A, ele_B, ele_C, i, ierr, nprocs, n_trisC, rank, serial_ele_A, serial_ele_B
  integer, parameter :: dim = 2, root = 0
  integer, dimension(:), allocatable :: ele_ownerA, ele_ownerB, unsA
  type(halo_type) :: halo
  type(intersections), dimension(:), allocatable :: map_AB
  type(tri_type) :: tri_A, tri_B
  type(tri_type), dimension(tri_buf_size) :: trisC
  type(vector_field) :: positionsA, positionsB
  real, parameter :: tol = 1.0e3 * epsilon(0.0)
  real :: t0, serial_time, parallel_time, serial_read_time, parallel_read_time
  logical :: fail

  real :: area_parallel, area_serial
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
  call write_parallel(buffer)

  ! Serial test
  if (rank == 0) then
    t0 = mpi_wtime()
    positionsA = read_triangle_files("data/square_0_02"//"_"//trim(nprocs_character), dim)
    positionsB = read_triangle_files("data/square_0_01"//"_"//trim(nprocs_character), dim)
    serial_ele_A = ele_count(positionsA)
    serial_ele_B = ele_count(positionsB)
    serial_read_time = mpi_wtime() - t0

    t0 = mpi_wtime()
    allocate(map_AB(serial_ele_A))
    call intersection_finder(positionsA%val, reshape(positionsA%mesh%ndglno, (/ele_loc(positionsA, 1), serial_ele_A/)), &
                         & positionsB%val, reshape(positionsB%mesh%ndglno, (/ele_loc(positionsB, 1), serial_ele_B/)), &
                         & map_AB)

    allocate(valsB(serial_ele_B))
    do ele_B = 1, serial_ele_B
      valsB(ele_B) = sum(positionsB%val(ele_nodes(positionsB, ele_B), 1)) / 3.0
    end do

    area_serial = 0.0
    do ele_A = 1, serial_ele_A
      tri_A%v = ele_val(positionsA, ele_A)

      do i = 1, map_AB(ele_A)%n
        ele_B = map_AB(ele_A)%v(i)
        tri_B%v = ele_val(positionsB, ele_B)

        call intersect_tris(tri_A, tri_B, trisC, n_trisC)

        do ele_C = 1, n_trisC
          area_serial = area_serial + valsB(ele_B) * triangle_area(trisC(ele_C)%v)
        end do
      end do
    end do

    deallocate(valsB)

    serial_time = mpi_wtime() - t0
    call deallocate(map_AB)
    deallocate(map_AB)

    call deallocate(positionsA)
    call deallocate(positionsB)
  end if

  call MPI_Barrier(MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)

  t0 = mpi_wtime()
  positionsA = read_triangle_files(trim("data/square_0_02_")//trim(nprocs_character)//"_"//trim(rank_character), dim)
  call read_halo("data/square_0_02"//"_"//trim(nprocs_character), halo, level = 2)
  allocate(ele_ownerA(ele_count(positionsA)))
  call element_ownership(node_count(positionsA), reshape(positionsA%mesh%ndglno, (/ele_loc(positionsA, 1), ele_count(positionsA)/)), halo, ele_ownerA)
  allocate(unsA(node_count(positionsA)))
  call universal_node_numbering(halo, unsA)
  call deallocate(halo)

  positionsB = read_triangle_files(trim("data/square_0_01_")//trim(nprocs_character)//"_"//trim(rank_character), dim)
  call read_halo("data/square_0_01"//"_"//trim(nprocs_character), halo, level = 2)
  allocate(ele_ownerB(ele_count(positionsB)))
  call element_ownership(node_count(positionsB), reshape(positionsB%mesh%ndglno, (/ele_loc(positionsB, 1), ele_count(positionsB)/)), halo, ele_ownerB)
  call deallocate(halo)
  parallel_read_time = mpi_wtime() - t0

  t0 = mpi_wtime()
  area_parallel = 0.0
  call parallel_supermesh(positionsA%val, reshape(positionsA%mesh%ndglno, (/ele_loc(positionsA, 1), ele_count(positionsA)/)), unsA, ele_ownerA, &
                      &   positionsB%val, reshape(positionsB%mesh%ndglno, (/ele_loc(positionsB, 1), ele_count(positionsB)/))      , ele_ownerB, &
                      &   local_donor_ele_data, local_intersection_calculation)
  parallel_time = mpi_wtime() - t0

  call MPI_Allreduce(area_parallel, area_parallel, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  call MPI_Allreduce(parallel_time, parallel_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  call MPI_Allreduce(parallel_read_time, parallel_read_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)

  call deallocate(positionsA)
  call deallocate(positionsB)
  deallocate(unsA)

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

    fail = fnequals(area_parallel, area_serial, tol = tol)
    call report_test("[test_parallel_partition_complete_ab areas]", fail, .FALSE., "Should give the same areas of intersection")
    if (fail) then
      print "(a,e25.17e3,a,e25.17e3,a)", ": Area, serial   =", area_serial, &
                                       & ": Area, parallel =", area_parallel
    end if
  end if

contains

  subroutine local_donor_ele_data(eles, data, ndata)
    use iso_c_binding, only : c_ptr
    implicit none
    integer, dimension(:), intent(in) :: eles
    type(c_ptr), intent(out)          :: data
    integer, intent(out)              :: ndata

  end subroutine local_donor_ele_data

  subroutine local_intersection_calculation(positions_c, ele_a, ele_b, data_b, ndata_b)
    use iso_c_binding, only : c_ptr
    implicit none
    ! dim x loc_c x nelements_c
    real, dimension(:, :, :), intent(in) :: positions_c
    integer, intent(in)     :: ele_a
    integer, intent(in)     :: ele_b
    type(c_ptr), intent(in) :: data_b
    integer, intent(in)     :: ndata_b

    integer :: ele_c

    do ele_c = 1, size(positions_c, 3)
      area_parallel = area_parallel + triangle_area(positions_c(:, :, ele_c))
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

end subroutine test_parallel_partition_complete_ab
