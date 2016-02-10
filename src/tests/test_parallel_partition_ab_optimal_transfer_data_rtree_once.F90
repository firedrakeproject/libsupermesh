#include "libsupermesh_debug.h"

subroutine test_parallel_partition_ab_optimal_transfer_data_rtree_once

  use iso_fortran_env, only : output_unit
  use iso_c_binding, only : c_double, c_int, c_int8_t, c_ptr, c_f_pointer, c_loc, c_size_t, c_sizeof

  use libsupermesh_unittest_tools
  use libsupermesh_supermesh
  use libsupermesh_fields
  use libsupermesh_read_triangle
  use libsupermesh_tri_intersection
  use libsupermesh_unittest_tools
  use libsupermesh_intersection_finder
  use libsupermesh_read_halos
  use libsupermesh_halo_ownership
  use libsupermesh_parallel_supermesh

  implicit none

  type pointer_real
    real, dimension(:), pointer :: p
  end type pointer_real

  type pointer_byte
    integer(kind = c_int8_t), dimension(:), pointer :: p
  end type pointer_byte

  type pointer_integer
    integer, dimension(:), pointer :: p
  end type pointer_integer

#include <finclude/petsc.h90>

  integer :: i, j, k, l, m, n, sends, recvs, ele_A, ele_B, ele_C, n_trisC, nprocs, &
       & rank, ierr, serial_ele_A, serial_ele_B, test_parallel_ele_A, &
       & test_parallel_ele_B, position, dp_extent, int_extent, icount
  integer :: local_sum_a, local_sum_b, triangles, &
       & serial_local_iter, serial_local_iter_actual, local_iter, local_iter_actual
  real :: t0, t1, area_serial, area_parallel, &
       & local_time, serial_read_time, serial_time, local_read_time, &
       & local_time_intersection_only, local_other_time, integral_serial, integral_parallel
  real, dimension(:), allocatable :: areas_parallel, integrals_parallel, times_parallel, times_intersection_only_parallel, other_times_parallel, times_read_time_parallel
  integer, dimension(:), allocatable :: iters_parallel, iter_actual_parallel

  integer, dimension(:), allocatable   :: number_of_elements_to_receive, number_of_elements_to_send
  integer, dimension(:,:), allocatable :: number_of_elements_and_nodes_to_receive, number_of_elements_and_nodes_to_send
  integer, dimension(:), allocatable   :: temp_elements_uns
  type(pointer_integer), dimension(:), allocatable  :: send_element_uns
  integer, dimension(:), allocatable   :: request_send, request_recv
  integer, dimension(:,:), allocatable :: status_send, status_recv
  type(pointer_byte), dimension(:), allocatable  :: recv_buffer, send_buffer
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
  integer, dimension(:), allocatable  :: test_parallel_ele_B_array

  integer, dimension(:), allocatable :: ele_ownerA, ele_ownerB, unsB
  type(halo_type) :: halo

  type(intersections), dimension(:), allocatable :: map_AB
  character(len = 255) :: hostname
  character(len = 2047) :: buffer

  integer :: nintersections
  integer, dimension(:, :), allocatable :: comm_enlist_B
  real, dimension(:, :), allocatable :: comm_coords_B

  integer, dimension(:), allocatable :: nodes
  integer, dimension(:,:), allocatable :: nodes_translation, nodes_translation_tmp
  type(pointer_integer), dimension(:), allocatable  :: send_nodes_conectivity, recv_nodes_conectivity
  type(pointer_real), dimension(:), allocatable  :: send_nodes_values, recv_nodes_values
  integer(kind = c_int8_t), dimension(:), pointer :: buffer_mpi
  integer                                :: buffer_size

  real, dimension(:), allocatable :: valsB

  ! Data
  integer(kind = c_int8_t), dimension(:), allocatable :: data
  integer, dimension(2) :: data_b_header
  real, dimension(:), allocatable :: data_b_data

  integer                         :: status(MPI_STATUS_SIZE)
  integer(kind = c_int8_t), dimension(:), allocatable :: ldata

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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Serial runtime test !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(rank == root) then
    t0 = mpi_wtime()
    positionsA = read_triangle_files("data/square_0_5"//"_"//trim(nprocs_character), dim)
    positionsB = read_triangle_files("data/square_0_9"//"_"//trim(nprocs_character), dim)
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
      valsB(ele_B) = sum(positionsB%val(1, ele_nodes(positionsB, ele_B))) / 3.0
    end do
!write(*,*) "serial size(valsB):",size(valsB),", valsB:",valsB
    area_serial = 0.0
    integral_serial = 0.0
    do ele_A = 1, serial_ele_A
      tri_A%v = ele_val(positionsA, ele_A)

      do i = 1, map_AB(ele_A)%n
        ele_B = map_AB(ele_A)%v(i)
        tri_B%v = ele_val(positionsB, ele_B)

        call intersect_tris(tri_A, tri_B, trisC, n_trisC)

        serial_local_iter = serial_local_iter + 1

        do ele_C = 1, n_trisC
          area_serial = area_serial + triangle_area(trisC(ele_C)%v)
          integral_serial = integral_serial + valsB(ele_B) * triangle_area(trisC(ele_C)%v)
          serial_local_iter_actual = serial_local_iter_actual + 1
        end do
      end do
    end do

    deallocate(valsB)
    call deallocate(map_AB)
    deallocate(map_AB)

    serial_time = mpi_wtime() - t0

    call deallocate(positionsA)
    call deallocate(positionsB)
  end if

  call MPI_Barrier(MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)

  allocate(test_parallel_ele_B_array(0:nprocs - 1))
  if(rank == root) then
    allocate(areas_parallel(0:nprocs - 1), integrals_parallel(0:nprocs - 1), times_parallel(0:nprocs - 1), other_times_parallel(0:nprocs - 1), times_read_time_parallel(0:nprocs - 1))
    allocate(iters_parallel(0:nprocs - 1))
    allocate(iter_actual_parallel(0:nprocs - 1))
    allocate(times_intersection_only_parallel(0:nprocs - 1))
  else
    allocate(areas_parallel(0), integrals_parallel(0), times_parallel(0), other_times_parallel(0), times_read_time_parallel(0))
    allocate(iters_parallel(0))
    allocate(iter_actual_parallel(0))
    allocate(times_intersection_only_parallel(0))
  end if

  area_parallel = 0.0
  integral_parallel = 0.0

  call MPI_Barrier(MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)

  t0 = mpi_wtime()
  positionsA = read_triangle_files(trim("data/square_0_5_")//trim(nprocs_character)//"_"//trim(rank_character), dim)
  call read_halo("data/square_0_5"//"_"//trim(nprocs_character), halo, level = 2)
  allocate(ele_ownerA(ele_count(positionsA)))
  call element_ownership(node_count(positionsA), reshape(positionsA%mesh%ndglno, (/ele_loc(positionsA, 1), ele_count(positionsA)/)), halo, ele_ownerA)
  test_parallel_ele_A = count(ele_ownerA == rank)
  call deallocate(halo)

  positionsB = read_triangle_files(trim("data/square_0_9_")//trim(nprocs_character)//"_"//trim(rank_character), dim)
  call read_halo("data/square_0_9"//"_"//trim(nprocs_character), halo, level = 2)
  allocate(ele_ownerB(ele_count(positionsB)))
  call element_ownership(node_count(positionsB), reshape(positionsB%mesh%ndglno, (/ele_loc(positionsB, 1), ele_count(positionsB)/)), halo, ele_ownerB)
  test_parallel_ele_B = count(ele_ownerB == rank)
  allocate(unsB(node_count(positionsB)))
  call universal_node_numbering(halo, unsB)
  call deallocate(halo)

  local_read_time = mpi_wtime() - t0

  t0 = mpi_wtime()
  allocate(valsB(test_parallel_ele_B))
  do ele_B = 1, test_parallel_ele_B
    valsB(ele_B) = sum(positionsB%val(1, ele_nodes(positionsB, ele_B))) / 3.0
  end do

  allocate(bbox_a(positionsA%dim,2))
  bbox_a = partition_bbox(positionsA, ele_ownerA, rank)
  allocate(bbox_b(positionsB%dim,2))
  bbox_b = partition_bbox(positionsB, ele_ownerB, rank)

  allocate(parallel_bbox_a(2, positionsA%dim, 0:nprocs-1))
  allocate(parallel_bbox_b(2, positionsB%dim, 0:nprocs-1))

  ! Bounding box all-to-all
  call MPI_Allgather(bbox_a, 4, MPI_DOUBLE_PRECISION, &
    & parallel_bbox_a, 4, MPI_DOUBLE_PRECISION,    &
    & MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)

  call MPI_Allgather(bbox_b, 4, MPI_DOUBLE_PRECISION, &
    & parallel_bbox_b, 4, MPI_DOUBLE_PRECISION,    &
    & MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)

  call MPI_Allgather(test_parallel_ele_B, 1, MPI_INTEGER, &
    & test_parallel_ele_B_array, 1, MPI_INTEGER,    &
    & MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)

  allocate(send_element_uns(0:nprocs-1))
  do i=0,size(send_element_uns(:))-1
    nullify(send_element_uns(i)%p)
  end do

  allocate(send_buffer(0:nprocs - 1))
  do i = 0, size(send_buffer(:)) - 1
    nullify(send_buffer(i)%p)
  end do

  allocate(recv_buffer(0:nprocs-1))
  do i=0,size(recv_buffer(:))-1
    nullify(recv_buffer(i)%p)
  end do

  allocate(send_nodes_conectivity(0:nprocs-1))
  do i=0,size(send_nodes_conectivity(:))-1
    nullify(send_nodes_conectivity(i)%p)
  end do

  allocate(recv_nodes_conectivity(0:nprocs-1))
  do i=0,size(recv_nodes_conectivity(:))-1
    nullify(recv_nodes_conectivity(i)%p)
  end do

  allocate(send_nodes_values(0:nprocs-1))
  do i=0,size(send_nodes_values(:))-1
    nullify(send_nodes_values(i)%p)
  end do

  allocate(recv_nodes_values(0:nprocs-1))
  do i=0,size(recv_nodes_values(:))-1
    nullify(recv_nodes_values(i)%p)
  end do

  call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, dp_extent, ierr); CHKERRQ(ierr)
  call MPI_TYPE_EXTENT(MPI_INTEGER, int_extent, ierr); CHKERRQ(ierr)

  allocate(number_of_elements_and_nodes_to_receive(2,0:nprocs - 1), number_of_elements_and_nodes_to_send(2,0:nprocs - 1))
  number_of_elements_and_nodes_to_receive = 0
  number_of_elements_and_nodes_to_send = 0

  allocate(partition_intersection_recv(0:nprocs - 1), partition_intersection_send(0:nprocs - 1))
  partition_intersection_recv = .false.
  partition_intersection_send = .false.
  sends = -1
  recvs = -1
  do i = 0, nprocs - 1
    l = 0

     if(partition_bboxes_intersect(bbox_a, parallel_bbox_b(:, :, i))) then
       if(i /= rank) then
         partition_intersection_recv(i) = .true.
       end if
     end if

    if(partition_bboxes_intersect(bbox_b, parallel_bbox_a(:, :, i))) then
      if(i /= rank) then
        partition_intersection_send(i) = .true.
        allocate(temp_elements_uns(ele_count(positionsB)))
        do ele_B = 1, ele_count(positionsB)
          if(ele_ownerB(ele_B) /= rank) cycle
          if(.not. partition_bboxes_intersect(bbox(ele_val(positionsB, ele_B)), parallel_bbox_a(:, :, i))) cycle
          l = l + 1                     ! Keep a counter
          temp_elements_uns(l) = ele_B  ! Keep the actual element
        end do

        allocate(nodes(l * 3))
        nodes = -12
        do k = 1, l
          nodes(k + (k-1)*positionsB%dim:k+((k-1)*positionsB%dim)+(positionsB%dim)) = ele_nodes(positionsB, temp_elements_uns(k))
        end do

        n = l/2
        allocate(nodes_translation(n,2))   ! May need to increase in size
        m = 0
        do j = minval(nodes), maxval(nodes)
          if (sum(nodes,mask=nodes.eq.j) == 0) cycle
          m = m + 1

          if (m < n) then
            nodes_translation(m,1) = m
            nodes_translation(m,2) = j
          else
            allocate(nodes_translation_tmp(m + n,2))
            n = m + n
            nodes_translation_tmp(1:size(nodes_translation,1),:) = nodes_translation
            deallocate(nodes_translation)
            allocate(nodes_translation(size(nodes_translation_tmp,1),2))
            nodes_translation = nodes_translation_tmp
            deallocate(nodes_translation_tmp)
            nodes_translation(m,1) = m
            nodes_translation(m,2) = j
          end if
        end do

        number_of_elements_and_nodes_to_send(1, i) = l
        number_of_elements_and_nodes_to_send(2, i) = m

!if (rank==2) write(*,*) rank,": i:",i,", nodes_translation(:,1):",nodes_translation(:,1),",nodes_translation(:,2):",nodes_translation(:,2)
!call flush()

        allocate(send_element_uns(i)%p(l))
        send_element_uns(i)%p = temp_elements_uns(1:l)

        allocate(send_nodes_conectivity(i)%p(number_of_elements_and_nodes_to_send(1, i) * (positionsB%dim + 1)))

        do k = 1, size(nodes)
          j = -1
          do n = 1, m
            if ( nodes(k) .eq. nodes_translation(n,2)) then
              j = nodes_translation(n,1)
              exit
            end if
          end do
          send_nodes_conectivity(i)%p(k) = j
        end do

        allocate(send_nodes_values(i)%p(number_of_elements_and_nodes_to_send(2,i) * (positionsB%dim)))
        send_nodes_values(i)%p = -22
        do n = 1, number_of_elements_and_nodes_to_send(2,i)
          send_nodes_values(i)%p(n + (n-1): n + (n-1) + 1) = node_val(positionsB, nodes_translation(n,2))
        end do

        ! ### Mesh send buffer allocation ###
! Removed IF, because sometimes we set partition_intersection_recv to TRUE;
! however, the following IF fails and we end up NOT sending anything.
! Now, we will send an extra message, so that the receiver (who already EXPECTS a message) can exit gracefully
        if ( (associated(send_element_uns(i)%p)) .and. (size(send_element_uns(i)%p) /= 0) ) then
          call local_donor_ele_data(send_element_uns(i)%p, data)

          buffer_size = int_extent + int_extent +                          &
                    &   (size(send_nodes_conectivity(i)%p) * int_extent) + &
                    &   (size(send_nodes_values(i)%p)      * dp_extent ) + &
                    &   int_extent + (size(data))

          allocate(buffer_mpi(buffer_size))
          position = 0

          call MPI_Pack(number_of_elements_and_nodes_to_send(1, i), &
               & 1,                                                 &
               & MPI_INTEGER, buffer_mpi, buffer_size, position, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
          call MPI_Pack(number_of_elements_and_nodes_to_send(2, i), &
               & 1,                                                 &
               & MPI_INTEGER, buffer_mpi, buffer_size, position, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
          call MPI_Pack(send_nodes_conectivity(i)%p,                &
               & size(send_nodes_conectivity(i)%p),                 &
               & MPI_INTEGER, buffer_mpi, buffer_size, position, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
          call MPI_Pack(send_nodes_values(i)%p,                     &
               & size(send_nodes_values(i)%p),                      &
               & MPI_DOUBLE_PRECISION, buffer_mpi, buffer_size, position, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
          call MPI_Pack(size(data),                                 &
               & 1,                                                 &
               & MPI_INTEGER, buffer_mpi, buffer_size, position, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
          call MPI_Pack(data,                                       &
               & size(data),                                        &
               & MPI_BYTE, buffer_mpi, buffer_size, position, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)

          deallocate(data)
        else
          buffer_size = int_extent + int_extent

          allocate(buffer_mpi(buffer_size))
          position = 0
          call MPI_Pack(number_of_elements_and_nodes_to_send(1, i), &
               & 1,                                                 &
               & MPI_INTEGER, buffer_mpi, buffer_size, position, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
          call MPI_Pack(number_of_elements_and_nodes_to_send(2, i), &
               & 1,                                                 &
               & MPI_INTEGER, buffer_mpi, buffer_size, position, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
        end if

        allocate(send_buffer(i)%p(buffer_size))
        send_buffer(i)%p = buffer_mpi
        sends = sends + 1

        call MPI_Isend(send_buffer(i)%p, &
               & size(send_buffer(i)%p), &
               & MPI_PACKED,             &
               & i, 0, MPI_COMM_WORLD, request_send(sends), ierr);  CHKERRQ(ierr)
        deallocate(buffer_mpi)

        deallocate(nodes, nodes_translation, temp_elements_uns)
      end if
    end if
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Parallel self-self runtime test !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  t1 = mpi_wtime()
  call rtree_intersection_finder_set_input(positionsA%val, reshape(positionsA%mesh%ndglno, (/ele_loc(positionsA, 1), ele_count(positionsA)/)))
  do ele_B = 1, ele_count(positionsB)
    if(ele_ownerB(ele_B) /= rank) cycle
    tri_B%v = ele_val(positionsB, ele_B)
    call rtree_intersection_finder_find(tri_B%v)
    call rtree_intersection_finder_query_output(nintersections)
    do i = 1, nintersections
      call rtree_intersection_finder_get_output(ele_A, i)
      if(ele_ownerA(ele_A) /= rank) cycle
      tri_A%v = ele_val(positionsA, ele_A)

      local_iter = local_iter + 1

      call intersect_tris(tri_A, tri_B, trisC, n_trisC)

      do ele_C = 1, n_trisC
        area_parallel = area_parallel + triangle_area(trisC(ele_C)%v)
        local_iter_actual = local_iter_actual + 1
        integral_parallel = integral_parallel + valsB(ele_B) * triangle_area(trisC(ele_C)%v)
      end do
    end do
  end do

  ! Receive the PACKED buffer and UNPACK
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Parallel self-other runtime test !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do i = 0,nprocs - 1
    if(rank == i) cycle

    if(partition_intersection_recv(i)) then
      call MPI_Probe(i, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr);  CHKERRQ(ierr)
      call MPI_Get_Count(status, MPI_PACKED, icount, ierr);  CHKERRQ(ierr)

      allocate(recv_buffer(i)%p(icount))
      call MPI_Recv (recv_buffer(i)%p, size(recv_buffer(i)%p), MPI_PACKED, & 
                &  i, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr);  CHKERRQ(ierr)

      buffer_size = size(recv_buffer(i)%p)
      position = 0
      call MPI_UnPack ( recv_buffer(i)%p, buffer_size,  &
           & position, m,                    &
           & 1,                              &
           & MPI_INTEGER, MPI_COMM_WORLD, IERR); CHKERRQ(ierr)
      allocate(recv_nodes_conectivity(i)%p(m * (positionsB%dim + 1)))

      call MPI_UnPack ( recv_buffer(i)%p, buffer_size,  &
           & position, n,                    &
           & 1,                              &
           & MPI_INTEGER, MPI_COMM_WORLD, IERR); CHKERRQ(ierr)
      allocate(recv_nodes_values(i)%p(n * (positionsB%dim)))
      number_of_elements_and_nodes_to_receive(1, i) = m
      number_of_elements_and_nodes_to_receive(2, i) = n

      if(number_of_elements_and_nodes_to_receive(1, i) <= 0) cycle

     call MPI_UnPack ( recv_buffer(i)%p, buffer_size,  &
          & position, recv_nodes_conectivity(i)%p,     &
          & size(recv_nodes_conectivity(i)%p),         &
          & MPI_INTEGER, MPI_COMM_WORLD, IERR); CHKERRQ(ierr)

     call MPI_UnPack ( recv_buffer(i)%p, buffer_size,  &
          & position, recv_nodes_values(i)%p,     &
          & size(recv_nodes_values(i)%p),         &
          & MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR); CHKERRQ(ierr)

      call MPI_UnPack ( recv_buffer(i)%p, buffer_size,  &
           & position, k,                    &
           & 1,                              &
           & MPI_INTEGER, MPI_COMM_WORLD, IERR); CHKERRQ(ierr)

      allocate(ldata(k))
      call MPI_UnPack ( recv_buffer(i)%p, buffer_size,  &
           & position, ldata,                    &
           & k,                              &
           & MPI_BYTE, MPI_COMM_WORLD, IERR); CHKERRQ(ierr)

      allocate(comm_coords_B(positionsB%dim,     number_of_elements_and_nodes_to_receive(2, i)), &
             & comm_enlist_B(positionsB%dim + 1, number_of_elements_and_nodes_to_receive(1, i)))

      call unpack_data_b(ldata)

      m = 1
      do j = 1, number_of_elements_and_nodes_to_receive(1, i)
        comm_enlist_B(:, j) = recv_nodes_conectivity(i)%p(j * (positionsB%dim + 1) - positionsB%dim: j * (positionsB%dim + 1) )
      end do
      m = 1
      do j = 1, number_of_elements_and_nodes_to_receive(2, i)
        do l = 1, positionsB%dim
          comm_coords_B(l, j) = recv_nodes_values(i)%p(m)
          m = m + 1
        end do
      end do

      do ele_B = 1, size(comm_enlist_B, 2)
        do j = 1, positionsB%dim + 1
          tri_B%v(:,j) = comm_coords_B(:, comm_enlist_B(j,ele_B))
        end do

        call rtree_intersection_finder_find(tri_B%v)
        call rtree_intersection_finder_query_output(nintersections)
        do k = 1, nintersections
          call rtree_intersection_finder_get_output(ele_A, k)
          if(ele_ownerA(ele_A) /= rank) cycle

          tri_A%v = ele_val(positionsA, ele_A)

          local_iter = local_iter + 1

          call intersect_tris(tri_A, tri_B, trisC, n_trisC)

          do ele_C = 1, n_trisC
            area_parallel = area_parallel + triangle_area(trisC(ele_C)%v)
            local_iter_actual = local_iter_actual + 1
            integral_parallel = integral_parallel + data_b_data(ele_b) * triangle_area(trisC(ele_C)%v)
          end do
        end do
      end do

      deallocate(comm_coords_B, &
               & comm_enlist_B)
      deallocate(ldata)
    end if
  end do
  call rtree_intersection_finder_reset()

  sends = sends + 1
  call MPI_Waitall(sends, request_send(0:sends), status_send(:,0:sends), ierr);  CHKERRQ(ierr)

  local_time = mpi_wtime() - t0                    ! total
  local_time_intersection_only = mpi_wtime() - t1  ! only intersection
  local_other_time = t1 - t0                       ! other
!write(output_unit, *) rank, ": area_parallel:",area_parallel
  ! Gather remote results:
  call MPI_Gather(area_parallel, 1, MPI_DOUBLE_PRECISION, &
    & areas_parallel, 1, MPI_DOUBLE_PRECISION, &
    & root, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  call MPI_Gather(integral_parallel, 1, MPI_DOUBLE_PRECISION, &
    & integrals_parallel, 1, MPI_DOUBLE_PRECISION, &
    & root, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  call MPI_Gather(local_time, 1, MPI_DOUBLE_PRECISION, &
    & times_parallel, 1, MPI_DOUBLE_PRECISION, &
    & root, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  call MPI_Gather(local_read_time, 1, MPI_DOUBLE_PRECISION, &
    & times_read_time_parallel, 1, MPI_DOUBLE_PRECISION, &
    & root, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  call MPI_Gather(local_other_time, 1, MPI_DOUBLE_PRECISION, &
    & other_times_parallel, 1, MPI_DOUBLE_PRECISION, &
    & root, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  call MPI_Gather(local_time_intersection_only, 1, MPI_DOUBLE_PRECISION, &
    & times_intersection_only_parallel, 1, MPI_DOUBLE_PRECISION, &
    & root, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  call MPI_Gather(local_iter, 1, MPI_INTEGER, &
    & iters_parallel, 1, MPI_INTEGER, &
    & root, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  call MPI_Gather(local_iter_actual, 1, MPI_INTEGER, &
    & iter_actual_parallel, 1, MPI_INTEGER, &
    & root, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  call MPI_Reduce(test_parallel_ele_A, local_sum_a, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  call MPI_Reduce(test_parallel_ele_B, local_sum_b, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)

  if(rank == root) then
    write(output_unit, "(i5,a,F19.15,a)") rank, ": Total serial read time          : ", serial_read_time," ."
    write(output_unit, "(i5,a,F19.15,a)") rank, ": Total serial intersection time  : ", serial_time," ."
    write(output_unit, "(a)") ""
    write(output_unit, "(i5,a,F19.15,a)") rank, ": Total parallel read time        : ", sum(times_read_time_parallel)," ."
    write(output_unit, "(i5,a,F19.15,a)") rank, ": Max parallel read time          : ", maxval(times_read_time_parallel)," ."
    write(output_unit, "(a)") ""
    write(output_unit, "(i5,a,F19.15,a)") rank, ": Total parallel intersection time: ", sum(times_intersection_only_parallel)," ."
    write(output_unit, "(i5,a,F19.15,a)") rank, ": Max parallel intersection time  : ", maxval(times_intersection_only_parallel)," ."
    write(output_unit, "(a)") ""
    write(output_unit, "(i5,a,F19.15,a)") rank, ": Total parallel time             : ", sum(times_parallel)," ."
    write(output_unit, "(i5,a,F19.15,a)") rank, ": Max parallel time               : ", maxval(times_parallel)," ."
    write(output_unit, "(a)") ""
    write(output_unit, "(i5,a,F19.15,a)") rank, ": Total other parallel time       : ", sum(other_times_parallel)," ."
    write(output_unit, "(i5,a,F19.15,a)") rank, ": Max other parallel time         : ", maxval(other_times_parallel)," ."
    write(output_unit, "(a)") ""
    write(output_unit, "(i5,a,F19.15,a)") rank, ": Total serial intersection area   : ", area_serial," ."
    write(output_unit, "(i5,a,F19.15,a)") rank, ": Total parallel intersection area : ", sum(areas_parallel)," ."
    write(output_unit, "(a)") ""
    write(output_unit, "(i5,a,F19.15,a)") rank, ": Total serial integral           : ", integral_serial," ."
    write(output_unit, "(i5,a,F19.15,a)") rank, ": Total parallel integral         : ", sum(integrals_parallel)," ."

    fail = fnequals(sum(areas_parallel), area_serial, tol = tol)
    call report_test("[test_parallel_partition_ab areas]", fail, .FALSE., "Should give the same areas of intersection")
    if (fail) then
      print "(a,e25.17e3,a,e25.17e3,a)", ": Total parallel intersection area:", sum(areas_parallel),&
                                       & ", total serial intersection area:",area_serial,"."
      write (*,*) "areas_parallel:",areas_parallel
      write (*,*) "iters_parallel:",iters_parallel
      print "(a,i15,a,i15,a)", "Serial iters:",serial_local_iter,", parallel iters:",sum(iters_parallel),"."
    end if

    fail = fnequals(sum(integrals_parallel), integral_serial, tol = tol)
    call report_test("[test_parallel_partition_ab integrals]", fail, .FALSE., "Should give the same values of integration")
    if (fail) then
      print "(a,e25.17e3,a,e25.17e3,a)", ": Total parallel integral :", sum(integrals_parallel),&
                                       & ", total serial integral   :",integral_serial,"."
      write (*,*) "integrals_parallel:",integrals_parallel
    end if    

    fail = ( sum(iter_actual_parallel) .ne. serial_local_iter_actual)
    call report_test("[test_parallel_partition_ab iterations]", fail, .FALSE., "Should give the same number of iterations")
    if (fail) then
      print "(a,i15,a,i15,a)", ": Total parallel actual iterations:", sum(iter_actual_parallel),&
                            & ", total serial actual iterations  :",serial_local_iter_actual,"."
    end if

    fail = ( serial_ele_A .ne. local_sum_a )
    call report_test("[test_parallel_partition_ab elementsA]", fail, .FALSE., "Should give the same number of elements for mesh A")
    if (fail) then
      print "(a,i5,a,i5,a)", ": Total parallel elements for mesh A:",local_sum_a,&
                           & ", total serial elements for mesh A:",serial_ele_A,"."
    end if

    fail = ( serial_ele_B .ne. local_sum_b )
    call report_test("[test_parallel_partition_ab elementsB]", fail, .FALSE., "Should give the same number of elements for mesh B")
    if (fail) then
      print "(a,i10,a,i10,a)", ": Total parallel elements for mesh B:",local_sum_b,&
                           & ", total serial elements for mesh B:",serial_ele_B,"."
    end if
  end if

  do i=0,size(send_element_uns(:))-1
    if ( .NOT. ASSOCIATED(send_element_uns(i)%p) ) cycle
    deallocate(send_element_uns(i)%p)
  end do
  deallocate(send_element_uns)

  do i=0,size(send_buffer(:))-1
    if ( .NOT. ASSOCIATED(send_buffer(i)%p) ) cycle
    deallocate(send_buffer(i)%p)
  end do
  deallocate(send_buffer)

  do i=0,size(recv_buffer(:))-1
    if ( .NOT. ASSOCIATED(recv_buffer(i)%p) ) cycle
    deallocate(recv_buffer(i)%p)
  end do
  deallocate(recv_buffer)

  do i=0,size(send_nodes_conectivity(:))-1
    if ( .NOT. ASSOCIATED(send_nodes_conectivity(i)%p) ) cycle
    deallocate(send_nodes_conectivity(i)%p)
  end do
  deallocate(send_nodes_conectivity)

  do i=0,size(recv_nodes_conectivity(:))-1
    if ( .NOT. ASSOCIATED(recv_nodes_conectivity(i)%p) ) cycle
    deallocate(recv_nodes_conectivity(i)%p)
  end do
  deallocate(recv_nodes_conectivity)

  do i=0,size(send_nodes_values(:))-1
    if ( .NOT. ASSOCIATED(send_nodes_values(i)%p) ) cycle
    deallocate(send_nodes_values(i)%p)
  end do
  deallocate(send_nodes_values)

  do i=0,size(recv_nodes_values(:))-1
    if ( .NOT. ASSOCIATED(recv_nodes_values(i)%p) ) cycle
    deallocate(recv_nodes_values(i)%p)
  end do
  deallocate(recv_nodes_values)

  deallocate(areas_parallel, integrals_parallel, times_parallel, times_read_time_parallel, other_times_parallel)
  deallocate(iters_parallel, iter_actual_parallel, times_intersection_only_parallel)
  deallocate(parallel_bbox_a, parallel_bbox_b, bbox_a, bbox_b, test_parallel_ele_B_array)
  deallocate(request_send, request_recv, status_send, status_recv, partition_intersection_send, partition_intersection_recv)
  call deallocate(positionsA)
  call deallocate(positionsB)
  deallocate(ele_ownerA, ele_ownerB)
  deallocate(unsB, valsB)
  deallocate(number_of_elements_to_receive, number_of_elements_to_send)

  call MPI_Barrier(MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)

contains

  subroutine local_donor_ele_data(eles, data)  
    integer, dimension(:), intent(in) :: eles
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

  subroutine unpack_data_b(data_b)
    integer(kind = c_int8_t), dimension(:), intent(in) :: data_b

    integer :: ierr, position

    call cleanup_data_b(data_b)

    position = 0
    call MPI_Unpack(data_b, size(data_b), position, data_b_header, 2, MPI_INTEGER, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    allocate(data_b_data(data_b_header(2)))
    call MPI_Unpack(data_b, size(data_b), position, data_b_data, size(data_b_data), MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)

  end subroutine unpack_data_b

  subroutine cleanup_data_b(data_b)
    integer(kind = c_int8_t), dimension(:), intent(in) :: data_b

    if(allocated(data_b_data)) deallocate(data_b_data)

  end subroutine cleanup_data_b

  subroutine intersection_calculation(positions_c, ele_a, ele_b, local)
    ! dim x loc_c x nelements_c
    real, dimension(:, :, :), intent(in) :: positions_c
    integer, intent(in)     :: ele_a
    integer, intent(in)     :: ele_b
    logical, intent(in)     :: local

    integer :: ele_c, nelements_c

    nelements_c = size(positions_c, 3)

!    do ele_c = 1, nelements_c
!      integral_parallel = integral_parallel + data_b_data(ele_b) * triangle_area(positions_c(:, :, ele_c))
!    end do

  end subroutine intersection_calculation

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

end subroutine test_parallel_partition_ab_optimal_transfer_data_rtree_once
