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
  use libsupermesh_parallel_supermesh, only : parallel_supermesh, finalise_parallel_supermesh

  implicit none

  type pointer_real
    real, dimension(:), pointer :: p
  end type pointer_real

  type pointer_integer
    integer, dimension(:), pointer :: p
  end type pointer_integer

#include <finclude/petsc.h90>

  integer :: i, j, k, l, sends, recvs, nnodes, ele_A, ele_B, ele_C, n_trisC, nprocs, &
       & rank, ierr, serial_ele_A, serial_ele_B, parallel_ele_A, &
       & parallel_ele_B
  integer :: local_sum_a, local_sum_b, triangles, &
       & serial_local_iter, serial_local_iter_actual, local_iter, local_iter_actual
  real :: t0, t1, area_serial, area_parallel, &
       & local_time, serial_read_time, serial_time, local_read_time, &
       & local_time_intersection_only, local_other_time
  real, dimension(:), allocatable :: areas_parallel, times_parallel, times_intersection_only_parallel, other_times_parallel, times_read_time_parallel
  integer, dimension(:), allocatable :: iters_parallel, iter_actual_parallel

  integer, dimension(:), allocatable   :: number_of_elements_to_receive, number_of_elements_to_send
  integer, dimension(:), allocatable   :: temp_elements_uns
  type(pointer_integer), dimension(:), allocatable  :: send_element_uns
  integer, dimension(:), allocatable   :: request_send, request_recv
  integer, dimension(:,:), allocatable :: status_send, status_recv
  type(pointer_real), dimension(:), allocatable  :: recv_buffer, send_buffer
  character(len = int(log10(real(huge(0)))) + 1) :: rank_character, nprocs_character
  type(tri_type) :: tri_A, tri_B
  type(tri_type), dimension(tri_buf_size) :: trisC
  type(vector_field) :: positionsA, positionsB
  integer, parameter :: dim = 2
  integer, parameter :: root = 0
  real, parameter :: tol = 1.0e3 * epsilon(0.0)
  logical :: fail
  logical, dimension(:), allocatable :: partition_intersection_recv, partition_intersection_send
  real, dimension(:,:), allocatable :: bbox_a, bbox_b
  real, dimension(:,:,:), allocatable :: parallel_bbox_a, parallel_bbox_b
  integer, dimension(:), allocatable  :: parallel_ele_B_array

  integer, dimension(:), allocatable :: ele_ownerA, ele_ownerB, unsA, unsB
  type(halo_type) :: halo

  type(intersections), dimension(:), allocatable :: map_AB
  character(len = 255) :: hostname
  character(len = 2047) :: buffer

! find out MY process ID, and how many processes were started.
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

  allocate(number_of_elements_to_receive(0:nprocs - 1), &
         & number_of_elements_to_send(0:nprocs -1),     &
         & request_send(0: nprocs - 1), request_recv(0: nprocs - 1), &
         & status_send(MPI_STATUS_SIZE, 0: nprocs - 1), status_recv(MPI_STATUS_SIZE, 0: nprocs - 1))
  number_of_elements_to_receive = -11
  number_of_elements_to_send    = -11
  forall(i = 0: nprocs - 1)
    status_send(:, i) = MPI_STATUS_IGNORE
    status_recv(:, i) = MPI_STATUS_IGNORE
  end forall
  request_send = MPI_REQUEST_NULL
  request_recv = MPI_REQUEST_NULL
  serial_local_iter = 0
  serial_local_iter_actual = 0
  local_iter = 0
  local_iter_actual = 0
  k = 0
  triangles = 0

  ! Serial test
  if (rank == 0) then
    t0 = mpi_wtime()
    positionsA = read_triangle_files("data/square_0_02"//"_"//trim(nprocs_character), dim)
    positionsB = read_triangle_files("data/square_0_01"//"_"//trim(nprocs_character), dim)
    serial_ele_A = ele_count(positionsA)
    serial_ele_B = ele_count(positionsB)
    serial_read_time = mpi_wtime() - t0

    t0 = mpi_wtime()
!   Use the intersection finder!!
    allocate(map_AB(serial_ele_A))
    call intersection_finder(positionsA%val, reshape(positionsA%mesh%ndglno, (/ele_loc(positionsA, 1), serial_ele_A/)), &
                         & positionsB%val, reshape(positionsB%mesh%ndglno, (/ele_loc(positionsB, 1), serial_ele_B/)), &
                         & map_AB)

    area_serial = 0.0
    do ele_A = 1, serial_ele_A
      tri_A%v = ele_val(positionsA, ele_A)

!      do ele_B = 1, serial_ele_B
      do i = 1, map_AB(ele_A)%n
        ele_B = map_AB(ele_A)%v(i)
      ! B. Use the libSuperMesh internal triangle intersector (using derived types as input)
        tri_B%v = ele_val(positionsB, ele_B)

        call intersect_tris(tri_A, tri_B, trisC, n_trisC)

        serial_local_iter = serial_local_iter + 1

        do ele_C = 1, n_trisC
          area_serial = area_serial + triangle_area(trisC(ele_C)%v)
          serial_local_iter_actual = serial_local_iter_actual + 1
        end do
      end do
    end do

    do i=1,size(map_AB)
      deallocate(map_AB(i)%v)
    end do
    deallocate(map_AB)

    serial_time = mpi_wtime() - t0

    call deallocate(positionsA)
    call deallocate(positionsB)
  end if

  call MPI_Barrier(MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)

  allocate(parallel_ele_B_array(0:nprocs - 1))
  if(rank == root) then
    allocate(areas_parallel(0:nprocs - 1), times_parallel(0:nprocs - 1), other_times_parallel(0:nprocs - 1), times_read_time_parallel(0:nprocs - 1))
    allocate(iters_parallel(0:nprocs - 1))
    allocate(iter_actual_parallel(0:nprocs - 1))
    allocate(times_intersection_only_parallel(0:nprocs - 1))
  else
    allocate(areas_parallel(0), times_parallel(0), other_times_parallel(0), times_read_time_parallel(0))
    allocate(iters_parallel(0))
    allocate(iter_actual_parallel(0))
    allocate(times_intersection_only_parallel(0))
  end if

  area_parallel = 0.0

  call MPI_Barrier(MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)

  t0 = mpi_wtime()
  positionsA = read_triangle_files(trim("data/square_0_02_")//trim(nprocs_character)//"_"//trim(rank_character), dim)
  call read_halo("data/square_0_02"//"_"//trim(nprocs_character), halo, level = 2)
  allocate(ele_ownerA(ele_count(positionsA)))
  call element_ownership(node_count(positionsA), reshape(positionsA%mesh%ndglno, (/ele_loc(positionsA, 1), ele_count(positionsA)/)), halo, ele_ownerA)
  parallel_ele_A = count(ele_ownerA == rank)
  allocate(unsA(node_count(positionsA)))
  call universal_node_numbering(halo, unsA)
  call deallocate(halo)

  positionsB = read_triangle_files(trim("data/square_0_01_")//trim(nprocs_character)//"_"//trim(rank_character), dim)
  call read_halo("data/square_0_01"//"_"//trim(nprocs_character), halo, level = 2)
  allocate(ele_ownerB(ele_count(positionsB)))
  call element_ownership(node_count(positionsB), reshape(positionsB%mesh%ndglno, (/ele_loc(positionsB, 1), ele_count(positionsB)/)), halo, ele_ownerB)
  parallel_ele_B = count(ele_ownerB == rank)
  allocate(unsB(node_count(positionsB)))
  call universal_node_numbering(halo, unsB)
  deallocate(unsB)
  call deallocate(halo)

  local_read_time = mpi_wtime() - t0

  t0 = mpi_wtime()
  call parallel_supermesh(positionsA%val, reshape(positionsA%mesh%ndglno, (/ele_loc(positionsA, 1), ele_count(positionsA)/)), unsA, ele_ownerA, &
                      &   positionsB%val, reshape(positionsB%mesh%ndglno, (/ele_loc(positionsB, 1), ele_count(positionsB)/))      , ele_ownerB, &
                      &   local_donor_ele_data, local_unpack_donor_ele_data, local_intersection_calculation)

  local_time = mpi_wtime() - t0                    ! total
  local_time_intersection_only = mpi_wtime() - t1  ! only intersection
  local_other_time = t1 - t0                       ! other


  call deallocate(positionsA)
  call deallocate(positionsB)

  call MPI_Barrier(MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  call finalise_parallel_supermesh()
  call cintersection_finder_reset(nnodes)

  if(rank == root) then
    write(output_unit, "(i5,a,F19.15,a)") rank, ": Total serial intersection area     : ", area_serial," ."
  end if

!   ! Gather remote results:
!   call MPI_Gather(area_parallel, 1, MPI_DOUBLE_PRECISION, &
!     & areas_parallel, 1, MPI_DOUBLE_PRECISION, &
!     & root, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
!   call MPI_Gather(local_time, 1, MPI_DOUBLE_PRECISION, &
!     & times_parallel, 1, MPI_DOUBLE_PRECISION, &
!     & root, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
!   call MPI_Gather(local_read_time, 1, MPI_DOUBLE_PRECISION, &
!     & times_read_time_parallel, 1, MPI_DOUBLE_PRECISION, &
!     & root, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
!   call MPI_Gather(local_other_time, 1, MPI_DOUBLE_PRECISION, &
!     & other_times_parallel, 1, MPI_DOUBLE_PRECISION, &
!     & root, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
!   call MPI_Gather(local_time_intersection_only, 1, MPI_DOUBLE_PRECISION, &
!     & times_intersection_only_parallel, 1, MPI_DOUBLE_PRECISION, &
!     & root, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
!   call MPI_Gather(local_iter, 1, MPI_INTEGER, &
!     & iters_parallel, 1, MPI_INTEGER, &
!     & root, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
!   call MPI_Gather(local_iter_actual, 1, MPI_INTEGER, &
!     & iter_actual_parallel, 1, MPI_INTEGER, &
!     & root, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
!   call MPI_Reduce(parallel_ele_A, local_sum_a, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
!   call MPI_Reduce(parallel_ele_B, local_sum_b, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
! 
!   if(rank == root) then
!     write(output_unit, "(i5,a,F19.15,a)") rank, ": Total serial read time          : ", serial_read_time," ."
!     write(output_unit, "(i5,a,F19.15,a)") rank, ": Total serial intersection time  : ", serial_time," ."
!     write(output_unit, "(a)") ""
!     write(output_unit, "(i5,a,F19.15,a)") rank, ": Total parallel read time        : ", sum(times_read_time_parallel)," ."
!     write(output_unit, "(i5,a,F19.15,a)") rank, ": Max parallel read time          : ", maxval(times_read_time_parallel)," ."
!     write(output_unit, "(a)") ""
!     write(output_unit, "(i5,a,F19.15,a)") rank, ": Total parallel intersection time: ", sum(times_intersection_only_parallel)," ."
!     write(output_unit, "(i5,a,F19.15,a)") rank, ": Max parallel intersection time  : ", maxval(times_intersection_only_parallel)," ."
!     write(output_unit, "(a)") ""
!     write(output_unit, "(i5,a,F19.15,a)") rank, ": Total parallel time             : ", sum(times_parallel)," ."
!     write(output_unit, "(i5,a,F19.15,a)") rank, ": Max parallel time               : ", maxval(times_parallel)," ."
!     write(output_unit, "(a)") ""
!     write(output_unit, "(i5,a,F19.15,a)") rank, ": Total other parallel time       : ", sum(other_times_parallel)," ."
!     write(output_unit, "(i5,a,F19.15,a)") rank, ": Max other parallel time         : ", maxval(other_times_parallel)," ."
!     write(output_unit, "(a)") ""
!     write(output_unit, "(i5,a,F19.15,a)") rank, ": Total serial intersection area     : ", area_serial," ."
!     write(output_unit, "(i5,a,F19.15,a)") rank, ": Total parallel intersection area   : ", sum(areas_parallel)," ."
! 
!     fail = fnequals(sum(areas_parallel), area_serial, tol = tol)
!     call report_test("[test_parallel_partition_ab areas]", fail, .FALSE., "Should give the same areas of intersection")
!     if (fail) then
!       print "(a,e25.17e3,a,e25.17e3,a)", ": Total parallel intersection area:", sum(areas_parallel),&
!                                        & ", total serial intersection area:",area_serial,"."
!       write (*,*) "areas_parallel:",areas_parallel
!       write (*,*) "iters_parallel:",iters_parallel
!       print "(a,i15,a,i15,a)", "Serial iters:",serial_local_iter,", parallel iters:",sum(iters_parallel),"."
!     end if
! 
!     fail = ( sum(iter_actual_parallel) .ne. serial_local_iter_actual)
!     call report_test("[test_parallel_partition_ab iterations]", fail, .FALSE., "Should give the same number of iterations")
!     if (fail) then
!       print "(a,i15,a,i15,a)", ": Total parallel actual iterations:", sum(iter_actual_parallel),&
!                             & ", total serial actual iterations  :",serial_local_iter_actual,"."
!     end if
! 
!     fail = ( serial_ele_A .ne. local_sum_a )
!     call report_test("[test_parallel_partition_ab elementsA]", fail, .FALSE., "Should give the same number of elements for mesh A")
!     if (fail) then
!       print "(a,i5,a,i5,a)", ": Total parallel elements for mesh A:",local_sum_a,&
!                            & ", total serial elements for mesh A:",serial_ele_A,"."
!     end if
! 
!     fail = ( serial_ele_B .ne. local_sum_b )
!     call report_test("[test_parallel_partition_ab elementsB]", fail, .FALSE., "Should give the same number of elements for mesh B")
!     if (fail) then
!       print "(a,i10,a,i10,a)", ": Total parallel elements for mesh B:",local_sum_b,&
!                            & ", total serial elements for mesh B:",serial_ele_B,"."
!     end if
!   end if

contains

  subroutine local_donor_ele_data(eles, data, ndata)
    use iso_c_binding, only : c_ptr
    implicit none
    integer, dimension(:), intent(in) :: eles
    type(c_ptr), intent(out)          :: data
    integer, intent(out)              :: ndata

  end subroutine local_donor_ele_data

  subroutine local_unpack_donor_ele_data(ele, proc, data, ndata, ele_data, nele_data)
    use iso_c_binding, only : c_ptr
    implicit none
    integer, intent(in)      :: ele
    integer, intent(in)      :: proc
    type(c_ptr), intent(in)  :: data
    integer, intent(in)      :: ndata
    type(c_ptr), intent(out) :: ele_data
    integer, intent(out)     :: nele_data

  end subroutine local_unpack_donor_ele_data

  subroutine local_intersection_calculation(positions_c, ele_a, proc_a, ele_b, ele_data_a, nele_data_a)
    use iso_c_binding, only : c_ptr
    implicit none
    ! dim x loc_c x nelements_c
    real, dimension(:, :, :), intent(in) :: positions_c
    integer, intent(in)     :: ele_a
    integer, intent(in)     :: proc_a
    integer, intent(in)     :: ele_b
    type(c_ptr), intent(in) :: ele_data_a
    integer, intent(in)     :: nele_data_a

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
