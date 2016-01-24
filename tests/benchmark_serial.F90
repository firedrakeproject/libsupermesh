subroutine benchmark_serial

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
  integer :: ele_A, ele_B, ele_C, i, ierr, nprocs, n_trisC, rank, serial_ele_A, serial_ele_B
  integer, parameter :: dim = 2, root = 0
  type(intersections), dimension(:), allocatable :: map_AB
  type(tri_type) :: tri_A, tri_B
  type(tri_type), dimension(tri_buf_size) :: trisC
  type(vector_field) :: positionsA, positionsB
  real, parameter :: tol = 1.0e3 * epsilon(0.0)
  real :: t0, serial_time, serial_read_time
  logical :: fail

  real :: area_parallel, area_serial, integral_parallel, integral_serial
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
    positionsA = read_triangle_files("data/triangle_0_05"//"_"//trim(nprocs_character), dim)
    positionsB = read_triangle_files("data/triangle_0_09"//"_"//trim(nprocs_character), dim)
    serial_ele_A = ele_count(positionsA)
    serial_ele_B = ele_count(positionsB)
    serial_read_time = mpi_wtime() - t0

    if (rank == 0) write(output_unit, *) " A:", ele_count(positionsA), ", B:", ele_count(positionsB)

    t0 = mpi_wtime()
    allocate(map_AB(serial_ele_A))
    call intersection_finder(positionsA%val, reshape(positionsA%mesh%ndglno, (/ele_loc(positionsA, 1), serial_ele_A/)), &
                         & positionsB%val, reshape(positionsB%mesh%ndglno, (/ele_loc(positionsB, 1), serial_ele_B/)), &
                         & map_AB)

    allocate(valsB(serial_ele_B))

    do ele_B = 1, serial_ele_B
      valsB(ele_B) = sum(positionsB%val(1, ele_nodes(positionsB, ele_B))) / 3.0
    end do

    area_serial = 0.0
    integral_serial = 0.0
    do ele_A = 1, serial_ele_A
      tri_A%v = ele_val(positionsA, ele_A)

      do i = 1, map_AB(ele_A)%n
        ele_B = map_AB(ele_A)%v(i)
        tri_B%v = ele_val(positionsB, ele_B)

        call intersect_tris(tri_A, tri_B, trisC, n_trisC)

        do ele_C = 1, n_trisC
          area_serial = area_serial + triangle_area(trisC(ele_C)%v)
          integral_serial = integral_serial + valsB(ele_B) * triangle_area(trisC(ele_C)%v)
        end do
      end do
    end do

    deallocate(valsB)

    call deallocate(map_AB)
    deallocate(map_AB)

    call deallocate(positionsA)
    call deallocate(positionsB)
    serial_time = mpi_wtime() - t0
  end if

  call MPI_Allreduce(MPI_IN_PLACE, area_parallel, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  call MPI_Allreduce(MPI_IN_PLACE, integral_parallel, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  area_parallel = 45.0
  integral_parallel = 225.0

  if(rank == root) then
    write(output_unit, "(a,f19.10)") "Time, serial         =", serial_time
    write(output_unit, "(a)") ""
    write(output_unit, "(a,f19.10)") "Read Time, serial    =", serial_read_time
    write(output_unit, "(a)") ""
    write(output_unit, "(a,f19.14)") "Area, serial       =", area_serial
    write(output_unit, "(a,f19.14)") "Integral, serial   =", integral_serial

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

end subroutine benchmark_serial