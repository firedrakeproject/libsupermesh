subroutine test_parallel_universal_node_numbering

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

  implicit none

  type pointer_real
    real, dimension(:), pointer :: p
  end type pointer_real

  type pointer_integer
    integer, dimension(:), pointer :: p
  end type pointer_integer

#include <finclude/petsc.h90>

  integer :: i, j, k, l, nprocs, rank, ierr, oldsize, &
       & serial_ele_A, serial_ele_B, parallel_ele_A, &
       & parallel_ele_B, meshA_nodes, meshB_nodes, no_of_same_nodes_A, no_of_same_nodes_B
  real :: t0, t1, area_parallel, serial_read_time
  real, dimension(:), allocatable :: times_parallel, other_times_parallel
  integer, dimension(:), allocatable :: meshA_nodes_parallel, meshB_nodes_parallel

  integer, dimension(:,:), allocatable :: unsA_parallel, unsB_parallel
  integer                              :: status(MPI_STATUS_SIZE)
  character(len = int(log10(real(huge(0)))) + 1) :: rank_character, nprocs_character
  type(vector_field) :: positionsA, positionsB
  integer, parameter :: dim = 2
  integer, parameter :: root = 0
  real, parameter :: tol = 1.0e3 * epsilon(0.0)
  logical :: fail
  real, dimension(:,:), allocatable :: meshA_node_val_parallel, meshB_node_val_parallel

  integer, dimension(:), allocatable :: ele_ownerA, ele_ownerB, unsA, unsB
  type(halo_type) :: halo

  character(len = 255)  :: hostname
  character(len = 2047) :: buffer

  INTEGER, PARAMETER    :: buffsize = 10000000
  CHARACTER             :: mpi_buffer ( buffsize )

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
  if ((rank == 0) .or. (mod(rank,10) == 0)) print *, trim(buffer)

  ! Serial test
  if (rank == 0) then
    t0 = mpi_wtime()
    positionsA = read_triangle_files("data/square_0_5"//"_"//trim(nprocs_character), dim)
    positionsB = read_triangle_files("data/square_0_09"//"_"//trim(nprocs_character), dim)
    serial_ele_A = ele_count(positionsA)
    serial_ele_B = ele_count(positionsB)
    call deallocate(positionsA)
    call deallocate(positionsB)
    serial_read_time = mpi_wtime() - t0
  end if

  call MPI_Barrier(MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)

  if(rank == root) then
    allocate(meshA_nodes_parallel(0:nprocs - 1), meshB_nodes_parallel(0:nprocs - 1), times_parallel(0:nprocs - 1), other_times_parallel(0:nprocs - 1))
  else
    allocate(meshA_nodes_parallel(0), meshB_nodes_parallel(0), times_parallel(0), other_times_parallel(0))
  end if

  meshA_nodes_parallel = 0
  meshB_nodes_parallel = 0

  area_parallel = 0.0

  call MPI_Buffer_Attach(mpi_buffer, buffsize, ierr);  CHKERRQ(ierr)
  call MPI_Barrier(MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)

  t0 = mpi_wtime()
  positionsA = read_triangle_files(trim("data/square_0_5_")//trim(nprocs_character)//"_"//trim(rank_character), dim)
  call read_halo("data/square_0_5"//"_"//trim(nprocs_character), halo, level = 2)
  allocate(ele_ownerA(ele_count(positionsA)))
  call element_ownership(node_count(positionsA), reshape(positionsA%mesh%ndglno, (/ele_loc(positionsA, 1), ele_count(positionsA)/)), halo, ele_ownerA)
  parallel_ele_A = count(ele_ownerA == rank)
  meshA_nodes = node_count(positionsA)
  allocate(unsA(node_count(positionsA)))
  call universal_node_numbering(halo, unsA)
  ! Gather remote results:
  call MPI_Gather(meshA_nodes, 1, MPI_INTEGER, &
    & meshA_nodes_parallel, 1, MPI_INTEGER, &
    & root, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)

  call MPI_Bsend(positionsA%val, meshA_nodes * 2, MPI_DOUBLE_PRECISION, 0, rank, &
          & MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  if (rank == 0) then
    allocate(meshA_node_val_parallel(maxval(meshA_nodes_parallel) * 2, 0:nprocs - 1))
    meshA_node_val_parallel = 0.0
    do i=0, nprocs - 1
      call MPI_Recv(meshA_node_val_parallel(:,i), meshA_nodes_parallel(i) * 2, MPI_DOUBLE_PRECISION, &
          & i, i, MPI_COMM_WORLD, status, ierr);  CHKERRQ(ierr)
    end do
  end if

  call MPI_Bsend(unsA, meshA_nodes, MPI_INTEGER, 0, rank, &
          & MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  if (rank == 0) then
    allocate(unsA_parallel(maxval(meshA_nodes_parallel), 0:nprocs - 1))
    unsA_parallel = 0
    do i=0, nprocs - 1
      call MPI_recv(unsA_parallel(:,i), meshA_nodes_parallel(i), MPI_INTEGER, &
          & i, i, MPI_COMM_WORLD, status, ierr);  CHKERRQ(ierr)
    end do
  end if
  call deallocate(halo)

  positionsB = read_triangle_files(trim("data/square_0_09_")//trim(nprocs_character)//"_"//trim(rank_character), dim)
  call read_halo("data/square_0_09"//"_"//trim(nprocs_character), halo, level = 2)
  allocate(ele_ownerB(ele_count(positionsB)))
  call element_ownership(node_count(positionsB), reshape(positionsB%mesh%ndglno, (/ele_loc(positionsB, 1), ele_count(positionsB)/)), halo, ele_ownerB)
  parallel_ele_B = count(ele_ownerB == rank)
  meshB_nodes = node_count(positionsB)
  allocate(unsB(node_count(positionsB)))
  call universal_node_numbering(halo, unsB)

  ! Gather remote results:
  call MPI_Gather(meshB_nodes, 1, MPI_INTEGER, &
    & meshB_nodes_parallel, 1, MPI_INTEGER, &
    & root, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)

  call MPI_Bsend(positionsB%val, meshB_nodes * 2, MPI_DOUBLE_PRECISION, 0, rank, &
          & MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  if (rank == 0) then
    allocate(meshB_node_val_parallel(maxval(meshB_nodes_parallel) * 2, 0:nprocs - 1))
    meshB_node_val_parallel = 0.0
    do i=0, nprocs - 1
      call MPI_Recv(meshB_node_val_parallel(:,i), meshB_nodes_parallel(i) * 2, MPI_DOUBLE_PRECISION, &
          & i, i, MPI_COMM_WORLD, status, ierr);  CHKERRQ(ierr)
    end do
  end if

  call MPI_Bsend(unsB, meshB_nodes, MPI_INTEGER, 0, rank, &
          & MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  if (rank == 0) then
    allocate(unsB_parallel(maxval(meshB_nodes_parallel), 0:nprocs - 1))
    unsB_parallel = 0
    do i=0, nprocs - 1
      call MPI_recv(unsB_parallel(:,i), meshB_nodes_parallel(i), MPI_INTEGER, &
          & i, i, MPI_COMM_WORLD, status, ierr);  CHKERRQ(ierr)
    end do
  end if

  if (rank == 0) then
    ! Check meshA
    no_of_same_nodes_A = 0
    fail = .false.
    do i=0, nprocs - 1
      do j=1, meshA_nodes_parallel(i)
        do k=0, nprocs - 1
          if (k == i) cycle
          do l=1, meshA_nodes_parallel(k)
            if (unsA_parallel(j,i) == unsA_parallel(l,k)) then
              fail = (fnequals((meshA_node_val_parallel(j*2 - 1,i)), (meshA_node_val_parallel(l*2 - 1,k)), tol = tol) .or. &
                   &  fnequals((meshA_node_val_parallel(j*2    ,i)), (meshA_node_val_parallel(l*2    ,k)), tol = tol)) .or. fail
              no_of_same_nodes_A = no_of_same_nodes_A + 1

              if (fail) then
                call report_test("[test_parallel_universal_node_numbering nodeA]", fail, .FALSE., "Should have the same co-ordinates.")
                print "(a,i5,a,i5,a,F19.15,a,F19.15,a,i5,a,i5,a,F19.15,a,F19.15)", &
                    & ": rank:",i,", node:",j,":",(meshA_node_val_parallel(j*2 - 1,i)), ",",(meshA_node_val_parallel(j*2,i)),&
                    & "  rank:",k,", node:",l,":",(meshA_node_val_parallel(l*2 - 1,k)), ",",(meshA_node_val_parallel(l*2,k))
              end if
            end if
          end do
        end do
      end do
    end do
    call report_test("[test_parallel_universal_node_numbering nodeA]", fail, .FALSE., "Should have the same co-ordinates.")

    ! Check meshB
    no_of_same_nodes_B = 0
    fail = .false.
    do i=0, nprocs - 1
      do j=1, meshB_nodes_parallel(i)
        do k=0, nprocs - 1
          if (k == i) cycle
          do l=1, meshB_nodes_parallel(k)
            if (unsB_parallel(j,i) == unsB_parallel(l,k)) then
              fail = (fnequals((meshB_node_val_parallel(j*2 - 1,i)), (meshB_node_val_parallel(l*2 - 1,k)), tol = tol) .or. &
                   &  fnequals((meshB_node_val_parallel(j*2    ,i)), (meshB_node_val_parallel(l*2    ,k)), tol = tol)) .or. fail
              no_of_same_nodes_B = no_of_same_nodes_B + 1

              if (fail) then
                call report_test("[test_parallel_universal_node_numbering nodeB]", fail, .FALSE., "Should have the same co-ordinates.")
                print "(a,i5,a,i5,a,F19.15,a,F19.15,a,i5,a,i5,a,F19.15,a,F19.15)", &
                    & ": rank:",i,", node:",j,":",(meshB_node_val_parallel(j*2 - 1,i)), ",",(meshB_node_val_parallel(j*2,i)),&
                    & "  rank:",k,", node:",l,":",(meshB_node_val_parallel(l*2 - 1,k)), ",",(meshB_node_val_parallel(l*2,k))
              end if
            end if
          end do
        end do
      end do
    end do
    call report_test("[test_parallel_universal_node_numbering nodeB]", fail, .FALSE., "Should have the same co-ordinates.")

  end if

  CALL MPI_Buffer_Detach ( mpi_buffer, oldsize, ierr);  CHKERRQ(ierr)
  call MPI_Barrier(MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)

  deallocate(unsA, unsB)
  call deallocate(halo)
  call deallocate(positionsA)
  call deallocate(positionsB)

contains

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

end subroutine test_parallel_universal_node_numbering
