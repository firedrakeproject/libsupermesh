subroutine test_parallel_partition_ab_optimal_transfer_data

  use iso_fortran_env, only : output_unit
  use iso_c_binding, only : c_double, c_int, c_int8_t, c_ptr, c_f_pointer, c_loc, c_size_t, c_sizeof

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

  type pointer_byte
    integer(kind = c_int8_t), dimension(:), pointer :: p
  end type pointer_byte

  type pointer_integer
    integer, dimension(:), pointer :: p
  end type pointer_integer

#include <finclude/petsc.h90>

  integer :: i, j, k, l, m, n, sends, recvs, nnodes, ele_A, ele_B, ele_C, n_trisC, nprocs, &
       & rank, ierr, serial_ele_A, serial_ele_B, test_parallel_ele_A, &
       & test_parallel_ele_B, position, dp_extent, int_extent
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

  integer :: nintersections, ntests
  integer, dimension(:, :), allocatable :: comm_enlist_B
  real, dimension(:, :), allocatable :: comm_coords_B

  integer, dimension(:), allocatable :: nodes, nodes_tmp
  integer, dimension(:,:), allocatable :: nodes_translation, nodes_translation_tmp
  type(pointer_integer), dimension(:), allocatable  :: send_nodes_conectivity, recv_nodes_conectivity
  type(pointer_real), dimension(:), allocatable  :: send_nodes_values, recv_nodes_values
  integer(kind = c_int8_t), dimension(:), pointer :: buffer_mpi
  integer                                :: buffer_size

  real, dimension(:), allocatable :: valsB

  ! Datatype(s)
  integer elements_nodes_datatype, oldtypes(0:1), blockcounts(0:1), offsets(0:1)
  
  ! Data
  integer(kind = c_int8_t), dimension(:), allocatable :: data
  integer, dimension(2) :: data_b_header
  real, dimension(:), allocatable :: data_b_data

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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Serial runtime test !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Serial test
  if(rank == root) then
    t0 = mpi_wtime()
    positionsA = read_triangle_files("data/square_0_2"//"_"//trim(nprocs_character), dim)
    positionsB = read_triangle_files("data/square_0_1"//"_"//trim(nprocs_character), dim)
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
  positionsA = read_triangle_files(trim("data/square_0_2_")//trim(nprocs_character)//"_"//trim(rank_character), dim)
  call read_halo("data/square_0_2"//"_"//trim(nprocs_character), halo, level = 2)
  allocate(ele_ownerA(ele_count(positionsA)))
  call element_ownership(node_count(positionsA), reshape(positionsA%mesh%ndglno, (/ele_loc(positionsA, 1), ele_count(positionsA)/)), halo, ele_ownerA)
  test_parallel_ele_A = count(ele_ownerA == rank)
  call deallocate(halo)

  positionsB = read_triangle_files(trim("data/square_0_1_")//trim(nprocs_character)//"_"//trim(rank_character), dim)
  call read_halo("data/square_0_1"//"_"//trim(nprocs_character), halo, level = 2)
  allocate(ele_ownerB(ele_count(positionsB)))
  call element_ownership(node_count(positionsB), reshape(positionsB%mesh%ndglno, (/ele_loc(positionsB, 1), ele_count(positionsB)/)), halo, ele_ownerB)
  test_parallel_ele_B = count(ele_ownerB == rank)
  allocate(unsB(node_count(positionsB)))
  call universal_node_numbering(halo, unsB)
!  if (rank == 1) write(output_unit, *) rank, ": size(unsb): ", size(unsB),", unsb: ", unsB
  call deallocate(halo)

  local_read_time = mpi_wtime() - t0

  t0 = mpi_wtime()
  allocate(valsB(test_parallel_ele_B))
  do i = 0, nprocs - 1
    call flush()
    call MPI_Barrier(MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  if (i == rank) then
write(*,*) ""
  do ele_B = 1, test_parallel_ele_B
    valsB(ele_B) = sum(positionsB%val(1, ele_nodes(positionsB, ele_B))) / 3.0
!if (i==2) write(*,*) rank,": ele_B:",ele_B,", unsB(ele_B):",unsB(ele_B),", valsB(ele_B):",valsB(ele_B)
  end do
!write(*,*) rank,": serial size(valsB):",size(valsB),", valsB:",valsB
  end if
    call flush()
    call MPI_Barrier(MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    call flush()
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

  call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, dp_extent, ierr)
  call MPI_TYPE_EXTENT(MPI_INTEGER, int_extent, ierr)

  ! Create Datatype
!  Setup description of the 2 MPI_INTEGER fields n, type 
   offsets(0) = 0
   oldtypes(0) = MPI_INTEGER
   blockcounts(0) = 2
!  Now define structured type and commit it 
  call MPI_TYPE_STRUCT(1, blockcounts, offsets, oldtypes,   &
                &      elements_nodes_datatype, ierr) ;  CHKERRQ(ierr)
  call MPI_TYPE_COMMIT(elements_nodes_datatype, ierr) ;  CHKERRQ(ierr)
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

    if(bboxes_intersect(bbox_a, parallel_bbox_b(:, :, i))) then
      if(i /= rank) then
        recvs = recvs + 1
        partition_intersection_recv(i) = .true.

!         call MPI_Irecv(number_of_elements_to_receive(i), 1, MPI_INTEGER, &
!             & i, i, MPI_COMM_WORLD, &
!             & request_recv(recvs), ierr);  CHKERRQ(ierr)
          call MPI_Irecv(number_of_elements_and_nodes_to_receive(:, i), 1, elements_nodes_datatype, &
              & i, i, MPI_COMM_WORLD, &
              & request_recv(recvs), ierr);  CHKERRQ(ierr)
      end if
    end if

    if(bboxes_intersect(bbox_b, parallel_bbox_a(:, :, i))) then
      if(i /= rank) then
        sends = sends + 1
        partition_intersection_send(i) = .true.
        allocate(temp_elements_uns(ele_count(positionsB)))
        do ele_B = 1, ele_count(positionsB)
          if(ele_ownerB(ele_B) /= rank) cycle
          if(.not. bboxes_intersect(bbox(ele_val(positionsB, ele_B)), parallel_bbox_a(:, :, i))) cycle
          l = l + 1                     ! Keep a counter
          temp_elements_uns(l) = ele_B  ! Keep the actual element
        end do

!if (rank == 1) write(output_unit, *) rank, ", l (elements): ", l," nodes (elements * 3):",l*3,", reals (nodes *2):",l*3*2
!if (rank == 1) write(output_unit, *) rank, ", l (elements): ", l," temp_elements_uns:",temp_elements_uns(1:l)
        allocate(nodes(l * 3))
        nodes = -12
        do k = 1, l
!if (rank == 1) write(output_unit, *) rank, ", k:",k,", element: ", temp_elements_uns(k)," nodes:",ele_nodes(positionsB, temp_elements_uns(k))
          nodes(k + (k-1)*2:((k-1)*2)+2) = ele_nodes(positionsB, temp_elements_uns(k))
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

!if (rank == 1) write(output_unit, *) rank, ": j:",j,", m (nodes):",m,", sum:",sum(nodes,mask=nodes.eq.j)/j,", translate:",nodes_translation(m,:)
        end do

!if (rank == 1) write(output_unit, *) rank, ", l (elements): ", l," nodes (            ):",m,", reals (nodes *2):",m*2
!if (rank == 1) write(output_unit, *) rank, ": nodes:",nodes
        number_of_elements_and_nodes_to_send(1, i) = l
        number_of_elements_and_nodes_to_send(2, i) = m
!         if (rank == 1) then
!           do j = 1, size(nodes)
!             n = -1
!             do k = 1, m
!               if (nodes_translation(k,2) .eq. nodes(j)) n = k
!             end do
!             write(output_unit, *) rank, ": old node:",nodes(j)," equal to new node:",nodes_translation(n,1)
!           end do
!         end if

!if (rank == 1) write(output_unit, *) rank, ": Sending to:",i,", number_of_elements_and_nodes_to_send:",number_of_elements_and_nodes_to_send(:, i)

        call MPI_Isend(number_of_elements_and_nodes_to_send(:, i), 1, elements_nodes_datatype, &
            & i, rank, MPI_COMM_WORLD, &
            & request_send(sends), ierr);  CHKERRQ(ierr)

!        number_of_elements_to_send(i) = l
        allocate(send_element_uns(i)%p(l))
        send_element_uns(i)%p = temp_elements_uns(1:l)
!if (rank == 1) write(output_unit, *) rank, ": Sending to:",i,", send_element_uns(i)%p:",send_element_uns(i)%p

        allocate(send_nodes_conectivity(i)%p(number_of_elements_and_nodes_to_send(1, i) * (positionsB%dim + 1)))
!if (rank == 1) write(output_unit, *) rank, ": Sending to:",i,", size of send_nodes_conectivity:",size(send_nodes_conectivity(i)%p)
        do k = 1, m
!if (rank == 1) write(output_unit, *) rank, ": k:",k,", translate:",nodes_translation(k,:)
        end do

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
!if ( (rank == 1) .and. (i==0)) write(output_unit, *) rank, ": Sending to:",i,", send_nodes_conectivity(i)%p:",send_nodes_conectivity(i)%p
!        send_nodes_conectivity(i)%p = -1 !nodes(1:12)!number_of_elements_and_nodes_to_send(2, i))
!        send_nodes_conectivity(i)%p = nodes(1:number_of_elements_and_nodes_to_send(2, i))


        allocate(send_nodes_values(i)%p(number_of_elements_and_nodes_to_send(2,i) * (positionsB%dim)))
        send_nodes_values(i)%p = -22
        do n = 1, number_of_elements_and_nodes_to_send(2,i)
!if (rank == 1) write(output_unit, *) rank, ": n:",n,", nodes_translation(n,2):",nodes_translation(n,2),", node_val:",node_val(positionsB, nodes_translation(n,2))
          send_nodes_values(i)%p(n + (n-1): n + (n-1) + 1) = node_val(positionsB, nodes_translation(n,2))
        end do
!if ( (rank == 1) .and. (i==0)) write(output_unit, *) rank, ": send_nodes_values(i)%p:",send_nodes_values(i)%p


!        do n = 1, number_of_elements_and_nodes_to_send(1,i)
!          tri_B%v = ele_val(positionsB, send_element_uns(i)%p(n))
!if (rank == 1) write(output_unit, *) rank, ": n:",n,", ele_nodes: ", ele_nodes(positionsB, send_element_uns(i)%p(n)),", (send_element_uns(j)%p(l):",send_element_uns(i)%p(n),", tri_B%v:",tri_B%v
!        end do

        deallocate(nodes)
        deallocate(nodes_translation)
        deallocate(temp_elements_uns)

        ! ### Mesh send buffer allocation ###
!        allocate(send_buffer(i)%p(l * (positionsB%dim + 1) * 2))
        allocate(send_buffer(i)%p((size(send_nodes_conectivity(i)%p) * int_extent) + &
                              &   (size(send_nodes_values(i)%p)      * dp_extent ) ) )
      end if
    end if

  end do

  recvs = recvs + 1
  call MPI_Waitall(recvs, request_recv(0:recvs), status_recv(:,0:recvs), ierr);  CHKERRQ(ierr)

  do i = 0, nprocs - 1
    if(partition_intersection_recv(i)) then
!if ((rank .eq. i) .or. (i .eq. 1)) write(output_unit, *) rank, ": i:",i,", number_of_elements_and_nodes_to_receive: ",  number_of_elements_and_nodes_to_receive(:, i)
    end if
  end do

!   call MPI_Barrier(MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
!   call flush()
!   call MPI_Barrier(MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)

  do i = 0, nprocs - 1
    if(partition_intersection_recv(i)) then
      ! ### Mesh receivebuffer allocation ###
      allocate(recv_buffer(i)%p((number_of_elements_and_nodes_to_receive(1, i) * (positionsB%dim + 1) * int_extent ) + &
                            &   (number_of_elements_and_nodes_to_receive(2, i) * (positionsB%dim) * dp_extent) ) )
      allocate(recv_nodes_conectivity(i)%p(number_of_elements_and_nodes_to_receive(1, i) * (positionsB%dim + 1)))
      allocate(recv_nodes_values(i)%p(number_of_elements_and_nodes_to_receive(2, i) * (positionsB%dim)))
    end if
  end do

!   call MPI_Barrier(MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
!   call flush()
!   call MPI_Barrier(MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
!   do j = 0, nprocs - 1
!     call MPI_Barrier(MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
!     call flush()
!     if (j == rank) then
!        do i = 0, nprocs - 1
!          if(partition_intersection_recv(i)) then
! write(output_unit, *) rank, ": i:",i,", elements:",size(send_nodes_conectivity(i)%p)/(positionsB%dim + 1),", nodes:",size(send_nodes_values(i)%p)/2,", send_buffer:", size(send_buffer(i)%p),", recv_buffer:",size(recv_buffer(i)%p)
! call flush()
!          end if
!       end do
!     end if
!     call flush()
!     call MPI_Barrier(MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
!     call flush()
!   end do

  j = -99
  do i = 0, nprocs - 1
    if(.not. associated(send_element_uns(i)%p)) cycle
    if(size(send_element_uns(i)%p) == 0) cycle

    call local_donor_ele_data(send_element_uns(i)%p, data)
    buffer_size = size(send_buffer(i)%p) + size(data)
    allocate(buffer_mpi(buffer_size))
    position = 0
    call MPI_Pack(send_nodes_conectivity(i)%p,      &
         & size(send_nodes_conectivity(i)%p),       &
         & MPI_INTEGER, buffer_mpi, buffer_size, position, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    call MPI_Pack(send_nodes_values(i)%p,           &
         & size(send_nodes_values(i)%p),            &
         & MPI_DOUBLE_PRECISION, buffer_mpi, buffer_size, position, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    call MPI_Pack(data,                             &
         & size(data),                              &
         & MPI_BYTE, buffer_mpi, buffer_size, position, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    deallocate(send_buffer(i)%p)
    send_buffer(i)%p => buffer_mpi

    deallocate(data)
  end do

  sends = sends + 1
  call MPI_Waitall(sends, request_send(0:sends), status_send(:,0:sends), ierr);  CHKERRQ(ierr)

  forall(i = 0: nprocs - 1)
    status_send(:, i) = MPI_STATUS_IGNORE
    status_recv(:, i) = MPI_STATUS_IGNORE
  end forall
  request_send = MPI_REQUEST_NULL
  request_recv = MPI_REQUEST_NULL

  sends = -1
  recvs = -1
  do i = 0,nprocs - 1
    if(rank == i) cycle

!    if(partition_intersection_recv(i) .and. (number_of_elements_to_receive(i) /=0) ) then
    if(partition_intersection_recv(i) .and. (number_of_elements_and_nodes_to_receive(2, i) /=0) ) then
      recvs = recvs + 1

!       call MPI_Irecv(recv_buffer(i)%p, &
!         & number_of_elements_to_receive(i) * (positionsB%dim + 1) * 2, &
!         & MPI_DOUBLE_PRECISION, &
!         & i, MPI_ANY_TAG, MPI_COMM_WORLD, &
!         & request_recv(recvs), ierr);  CHKERRQ(ierr)
      call MPI_Irecv(recv_buffer(i)%p, &
        & size(recv_buffer(i)%p), &
        & MPI_PACKED, &
        & i, MPI_ANY_TAG, MPI_COMM_WORLD, &
        & request_recv(recvs), ierr);  CHKERRQ(ierr)
    end if

    if( (partition_intersection_send(i)) .and. (size(send_buffer(i)%p) /=0) ) then
      sends = sends + 1

!       call MPI_Isend(send_buffer(i)%p, &
!         & size(send_buffer(i)%p), &
!         & MPI_DOUBLE_PRECISION, &
!         & i, 0, MPI_COMM_WORLD, request_send(sends), ierr);  CHKERRQ(ierr)
      call MPI_Isend(send_buffer(i)%p, &
        & size(send_buffer(i)%p), &
        & MPI_PACKED, &
        & i, 0, MPI_COMM_WORLD, request_send(sends), ierr);  CHKERRQ(ierr)
    end if
  end do

!write(output_unit, *) rank, ": positionsB%val:",positionsB%val,", positionsB%mesh%ndglno:",positionsB%mesh%ndglno

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Parallel self-self runtime test !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  t1 = mpi_wtime()
  call rtree_intersection_finder_set_input(positionsB%val, reshape(positionsB%mesh%ndglno, (/ele_loc(positionsB, 1), ele_count(positionsB)/)))
  do ele_A = 1, ele_count(positionsA)
    if(ele_ownerA(ele_A) /= rank) cycle
    tri_A%v = ele_val(positionsA, ele_A)
    call rtree_intersection_finder_find(tri_A%v)
    call rtree_intersection_finder_query_output(nintersections)
    do i = 1, nintersections
      call rtree_intersection_finder_get_output(ele_B, i)
      if(ele_ownerB(ele_B) /= rank) cycle
      tri_B%v = ele_val(positionsB, ele_B)

      local_iter = local_iter + 1

      call intersect_tris(tri_A, tri_B, trisC, n_trisC)

      do ele_C = 1, n_trisC
        area_parallel = area_parallel + triangle_area(trisC(ele_C)%v)
        local_iter_actual = local_iter_actual + 1
      end do
    end do
  end do
  call rtree_intersection_finder_reset(ntests)

  sends = sends + 1
  recvs = recvs + 1
  call MPI_Waitall(sends, request_send(0:sends), status_send(:,0:sends), ierr);  CHKERRQ(ierr)
  call MPI_Waitall(recvs, request_recv(0:recvs), status_recv(:,0:recvs), ierr);  CHKERRQ(ierr)


  do j = 0, nprocs - 1
!     call MPI_Barrier(MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
!     call flush()
    if (j == rank) then
       do i = 0, nprocs - 1
         if(partition_intersection_recv(i)) then

     position = 0
     buffer_size = size(recv_buffer(i)%p)
     call MPI_UnPack ( recv_buffer(i)%p, buffer_size,  &
          & position, recv_nodes_conectivity(i)%p,     &
          & size(recv_nodes_conectivity(i)%p),         &
          & MPI_INTEGER, MPI_COMM_WORLD, IERR); CHKERRQ(ierr)
     call MPI_UnPack ( recv_buffer(i)%p, buffer_size,  &
          & position, recv_nodes_values(i)%p,     &
          & size(recv_nodes_values(i)%p),         &
          & MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR); CHKERRQ(ierr)

!write(output_unit, *) ""
!write(output_unit, *) rank, ": i:",i,", number_of_elements_to_receive(i): ", number_of_elements_to_receive(i)
!write(output_unit, *) rank, ": i:",i,", (ELEMENTS) number_of_elements_and_nodes_to_receive(1, i): ", number_of_elements_and_nodes_to_receive(1, i),", (NODES) number_of_elements_and_nodes_to_receive(2, i): ", number_of_elements_and_nodes_to_receive(2, i)
!write(output_unit, *) rank, ": i:",i,", size(recv_nodes_conectivity(i)%p): ", size(recv_nodes_conectivity(i)%p),", size(recv_nodes_values(i)%p): ", size(recv_nodes_values(i)%p)

!write(output_unit, *) rank, ": i:",i,", size(send_nodes_conectivity(i)%p): ", size(send_nodes_conectivity(i)%p),", send_nodes_conectivity(i)%p:",send_nodes_conectivity(i)%p
!if ( (rank == 0) .and. (i==1)) write(output_unit, *) rank, ": i:",i,", size(recv_nodes_conectivity(i)%p): ", size(recv_nodes_conectivity(i)%p),", recv_nodes_conectivity(i)%p:",recv_nodes_conectivity(i)%p
! write(output_unit, *) rank, ": i:",i,", size(send_nodes_values(i)%p): ", size(send_nodes_values(i)%p),", send_nodes_values(i)%p:",send_nodes_values(i)%p
!if ( (rank == 0) .and. (i==1))  write(output_unit, *) rank, ": i:",i,", size(recv_nodes_values(i)%p): ", size(recv_nodes_values(i)%p),", recv_nodes_values(i)%p:",recv_nodes_values(i)%p
         end if
      end do
    end if
!     call flush()
!     call MPI_Barrier(MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
!     call flush()
  end do

!   call flush()
!   call MPI_Barrier(MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
!   call flush()

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Parallel self-other runtime test !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do i = 0, nprocs - 1
    if(i == rank) cycle
    if(number_of_elements_and_nodes_to_receive(1, i) <= 0) cycle
    allocate(comm_coords_B(positionsB%dim,     number_of_elements_and_nodes_to_receive(2, i)), &
           & comm_enlist_B(positionsB%dim + 1, number_of_elements_and_nodes_to_receive(1, i)))
!     m = 1
!     do j = 1, number_of_elements_to_receive(1, i)
!       do k = 1, positionsB%dim + 1
!         do l = 1, positionsB%dim
!           comm_coords_B(l, (j - 1) * (positionsB%dim + 1) + k) = recv_buffer(i)%p(positionsB%dim * (m - 1) + l)
!         end do
!         comm_enlist_B(k, j) = m
!         m = m + 1
!       end do
!     end do
!    call rtree_intersection_finder_set_input(comm_coords_B, comm_enlist_B)
!    call rtree_intersection_finder_set_input(recv_nodes_values(i)%p, recv_nodes_conectivity(i)%p)

!    do n = 1, number_of_elements_and_nodes_to_send(2,i)
!        send_nodes_values(i)%p(n + (n-1): n + (n-1) + 1) = node_val(positionsB, nodes_translation(n,2))
      
!    end do
    m = 1
    do j = 1, number_of_elements_and_nodes_to_receive(1, i)
!       do k = 1, positionsB%dim + 1
!         do l = 1, positionsB%dim
! !          comm_coords_B(l,:) = recv_nodes_values(i)%p()
! !          comm_coords_B(l, (j - 1) * (positionsB%dim + 1) + k) = recv_nodes_values(i)%p(j + (j-1))
!         end do
!         m = m + 1
!       end do
      comm_enlist_B(:, j) = recv_nodes_conectivity(i)%p(j * (positionsB%dim + 1) - positionsB%dim: j * (positionsB%dim + 1) )
    end do
    m = 1
    do j = 1, number_of_elements_and_nodes_to_receive(2, i)
      do l = 1, positionsB%dim
!      send_nodes_values(i)%p(n + (n-1): n + (n-1) + 1)
        comm_coords_B(l, j) = recv_nodes_values(i)%p(m)
        m = m + 1
      end do
    end do

!if ( (rank == 0) .and. (i==1))  write(output_unit, *) rank, ": i:",i,", size(comm_enlist_B): ", size(comm_enlist_B),", comm_enlist_B:",comm_enlist_B

!if ( (rank == 0) .and. (i==1))  write(output_unit, *) rank, ": i:",i,", size(comm_coords_B): ", size(comm_coords_B),", comm_coords_B:",comm_coords_B

   call rtree_intersection_finder_set_input(comm_coords_B, comm_enlist_B)
    do ele_A = 1, ele_count(positionsA)
      if(ele_ownerA(ele_A) /= rank) cycle
      tri_A%v = ele_val(positionsA, ele_A)
      call rtree_intersection_finder_find(tri_A%v)
      call rtree_intersection_finder_query_output(nintersections)
      do k = 1, nintersections
        call rtree_intersection_finder_get_output(ele_B, k)
!if ( (rank == 0) .and. (i==1))  write(output_unit, *) rank, ": i:",i,", ele_B:",ele_B
!call flush()
!call flush()
!if ( (rank == 0) .and. (i==1))  write(output_unit, *) rank, ": comm_enlist_B(ele_B):",comm_enlist_B(:,ele_B)
!call flush()
!call flush()
        do j = 1, positionsB%dim + 1
!          tri_B%v = comm_coords_B(:, (ele_B - 1) * (positionsB%dim + 1) + 1:ele_B * (positionsB%dim + 1))
          tri_B%v(:,j) = comm_coords_B(:, comm_enlist_B(j,ele_B))
!if ( (rank == 0) .and. (i==1))  write(output_unit, *) rank, ": comm_coords_B(:, comm_enlist_B(j,ele_B) * (positionsB%dim + 1) + 1:comm_enlist_B(j,ele_B) * (positionsB%dim + 1)):",comm_coords_B(:, comm_enlist_B(j,ele_B))
        end do
!if ( (rank == 0) .and. (i==1))  write(output_unit, *) rank, ": k:",k,", tri_B%v: ", tri_B%v
!call flush()

        local_iter = local_iter + 1

        call intersect_tris(tri_A, tri_B, trisC, n_trisC)

        do ele_C = 1, n_trisC
          area_parallel = area_parallel + triangle_area(trisC(ele_C)%v)
          local_iter_actual = local_iter_actual + 1
        end do
      end do
    end do
    call rtree_intersection_finder_reset(ntests)
    deallocate(comm_coords_B, &
             & comm_enlist_B)
 end do
  local_time = mpi_wtime() - t0                    ! total
  local_time_intersection_only = mpi_wtime() - t1  ! only intersection
  local_other_time = t1 - t0                       ! other
write(output_unit, *) rank, ": area_parallel:",area_parallel
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
    write(output_unit, "(i5,a,F19.15,a)") rank, ": Total serial intersection area     : ", area_serial," ."
    write(output_unit, "(i5,a,F19.15,a)") rank, ": Total parallel intersection area   : ", sum(areas_parallel)," ."
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

  call MPI_TYPE_FREE(elements_nodes_datatype, ierr);  CHKERRQ(ierr)

  deallocate(areas_parallel, integrals_parallel, times_parallel, times_read_time_parallel, other_times_parallel)
  deallocate(iters_parallel, iter_actual_parallel, times_intersection_only_parallel)
  deallocate(parallel_bbox_a, parallel_bbox_b, bbox_a, bbox_b, test_parallel_ele_B_array)
  deallocate(request_send, request_recv, status_send, status_recv, partition_intersection_send, partition_intersection_recv)
  call deallocate(positionsA)
  call deallocate(positionsB)
  deallocate(ele_ownerA, ele_ownerB)
  deallocate(unsB, valsB)
  deallocate(number_of_elements_to_receive, number_of_elements_to_send)

  call cintersection_finder_reset(nnodes)

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

  subroutine intersection_calculation(positions_c, ele_a, ele_b)
    ! dim x loc_c x nelements_c
    real, dimension(:, :, :), intent(in) :: positions_c
    integer, intent(in)     :: ele_a
    integer, intent(in)     :: ele_b
    
    integer :: ele_c, nelements_c
    
    nelements_c = size(positions_c, 3)
    
    do ele_c = 1, nelements_c
      integral_parallel = integral_parallel + data_b_data(ele_b) * triangle_area(positions_c(:, :, ele_c))
    end do
    
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

end subroutine test_parallel_partition_ab_optimal_transfer_data
