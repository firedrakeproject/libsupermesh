#include "libsupermesh_debug.h"

module libsupermesh_parallel_supermesh

  use iso_c_binding, only : c_int8_t
  
  use libsupermesh_debug
  use libsupermesh_debug_parameters, only : debug_log_unit
  use libsupermesh_integer_hash_table
  use libsupermesh_read_halos
  use libsupermesh_intersection_finder
  use libsupermesh_tri_intersection
  use libsupermesh_tet_intersection
  use libsupermesh_supermesh
  use libsupermesh_integer_set

  implicit none

#include <mpif.h>

  private

  public :: parallel_supermesh, print_times
  public :: bbox, partition_bbox, partition_bboxes_intersect

  type buffer
    integer(kind = c_int8_t), dimension(:), pointer :: v => null()
  end type buffer
  
  interface partition_bbox
    module procedure partition_bbox_real, partition_bbox_vector_field
  end interface partition_bbox

  logical, save :: parallel_supermesh_allocated = .false.
  
  real, dimension(:, :), allocatable, save :: bbox_a, bbox_b
  real, dimension(:, :, :), allocatable, save :: parallel_bbox_a, parallel_bbox_b

  integer :: nsends
  integer, dimension(:), allocatable, save :: send_requests
  logical, dimension(:), allocatable, save :: recv
  type(buffer), dimension(:), allocatable, save :: send_buffer

  integer, save :: mpi_comm, nprocs, rank
  real :: t0, all_to_all, point_to_point, local_compute, remote_compute, t1

contains

  subroutine parallel_supermesh(positions_a, enlist_a, ele_owner_a, &
                             &  positions_b, enlist_b, ele_owner_b, &
                             &  donor_ele_data, unpack_data_b, intersection_calculation, &
                             &  comm)
    ! dim x nnodes_a
    real, dimension(:, :), intent(in)    :: positions_a
    ! loc_a x nelements_a
    integer, dimension(:, :), intent(in) :: enlist_a
    ! elements_a
    integer, dimension(:), intent(in)    :: ele_owner_a
    ! dim x nnodes_b
    real, dimension(:, :), intent(in)    :: positions_b
    ! loc_b x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    ! elements_b
    integer, dimension(:), intent(in)    :: ele_owner_b
    integer, optional, intent(in) :: comm
    interface
      subroutine donor_ele_data(nodes, eles, data)
        use iso_c_binding, only : c_int8_t
        implicit none
        integer, dimension(:), intent(in)                                :: nodes
        integer, dimension(:), intent(in)                                :: eles
        integer(kind = c_int8_t), dimension(:), allocatable, intent(out) :: data
      end subroutine donor_ele_data

      subroutine unpack_data_b(data_b)
        use iso_c_binding, only : c_int8_t
        implicit none
        integer(kind = c_int8_t), dimension(:), intent(in) :: data_b
      end subroutine unpack_data_b
      
      subroutine intersection_calculation(positions_a, positions_b, positions_c, nodes_b, ele_a, ele_b, local)
        use iso_c_binding, only : c_int8_t
        implicit none
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
      end subroutine intersection_calculation
    end interface

    ! 0. Allocate and initialise arrays
    call initialise_parallel_supermesh(comm)

    ! 1. Calculate and communicate bounding box data for donor and target
    call step_1(positions_a, enlist_a, ele_owner_a, positions_b, enlist_b, ele_owner_b)

    ! 2. Use bounding box data to cull donor mesh
    call step_2(positions_b, enlist_b, ele_owner_b, donor_ele_data)

    ! 5. Supermesh and call user specified element unpack and calculation functions
    call step_3(positions_a, enlist_a, ele_owner_a, positions_b, enlist_b, ele_owner_b, &
             &  unpack_data_b, intersection_calculation)

    ! 6. Deallocate
    call finalise_parallel_supermesh()

  end subroutine parallel_supermesh

  subroutine step_1(positions_a, enlist_a, ele_owner_a, positions_b, enlist_b, ele_owner_b)
    ! dim x nnodes_a
    real, dimension(:, :), intent(in) :: positions_a
    ! loc_a x nelements_a
    integer, dimension(:, :), intent(in) :: enlist_a
    ! elements_a
    integer, dimension(:), intent(in) :: ele_owner_a
    ! dim x nnodes_b
    real, dimension(:, :), intent(in) :: positions_b
    ! loc_b x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    ! elements_b
    integer, dimension(:), intent(in) :: ele_owner_b

    integer :: dim, ierr

    dim = size(positions_a, 1)
    
    allocate(bbox_a(2, dim))
    bbox_a = partition_bbox(positions_a, enlist_a, ele_owner_a, rank)

    allocate(bbox_b(2, dim))
    bbox_b = partition_bbox(positions_b, enlist_b, ele_owner_b, rank)

    allocate(parallel_bbox_a(2, dim, nprocs))
    allocate(parallel_bbox_b(2, dim, nprocs))

#ifdef PROFILE
    t0 = mpi_wtime()
#else
    point_to_point = 0.0D0
    all_to_all = 0.0D0
    t0 = 0.0D0
    local_compute = 0.0D0
    remote_compute = 0.0D0
    t1 = 0.0D0
#endif
    call MPI_Allgather(bbox_a, 2 * dim, MPI_DOUBLE_PRECISION, &
      & parallel_bbox_a, 2 * dim, MPI_DOUBLE_PRECISION, &
      & mpi_comm, ierr);  assert(ierr == MPI_SUCCESS)

    call MPI_Allgather(bbox_b, 2 * dim, MPI_DOUBLE_PRECISION, &
      & parallel_bbox_b, 2 * dim, MPI_DOUBLE_PRECISION, &
      & mpi_comm, ierr);  assert(ierr == MPI_SUCCESS)
#ifdef PROFILE
    all_to_all = mpi_wtime() - t0
#endif

  end subroutine step_1


  subroutine step_2(positions_b, enlist_b, ele_owner_b, donor_ele_data)
    ! dim x nnodes_b
    real, dimension(:, :), intent(in)    :: positions_b
    ! loc_b x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    ! elements_b
    integer, dimension(:), intent(in)    :: ele_owner_b

    integer    :: ele_B, i, j, k, ierr, buffer_size, dp_extent, int_extent, position
    
    integer, dimension(:), allocatable :: elements_send, send_nodes_array
    integer :: nelements_send
    type(integer_set) :: nodes_send
    type(integer_hash_table) :: node_map
    integer :: nnodes_send
    integer :: node

    interface
      subroutine donor_ele_data(nodes, eles, data)
        use iso_c_binding, only : c_int8_t
        implicit none
        integer, dimension(:), intent(in)                                :: nodes
        integer, dimension(:), intent(in)                                :: eles
        integer(kind = c_int8_t), dimension(:), allocatable, intent(out) :: data
      end subroutine donor_ele_data
    end interface

    integer, dimension(:), allocatable :: send_enlist
    integer(kind = c_int8_t), dimension(:), allocatable :: data
    real, dimension(:), allocatable :: send_positions

#ifdef PROFILE
    t0 = -1
#endif

    call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, dp_extent, ierr)
    if(ierr /= MPI_SUCCESS) then
      libsupermesh_abort("MPI_TYPE_EXTENT double precision error")
    end if
    call MPI_TYPE_EXTENT(MPI_INTEGER, int_extent, ierr)
    if(ierr /= MPI_SUCCESS) then
      libsupermesh_abort("MPI_TYPE_EXTENT integer error")
    end if

    nsends = 0
    do i = 0, nprocs - 1
      if(i == rank) then
        recv(i + 1) = .false.
        cycle
      end if

      recv(i + 1) = partition_bboxes_intersect(bbox_a, parallel_bbox_b(:, :, i + 1))

      if(partition_bboxes_intersect(bbox_b, parallel_bbox_a(:, :, i + 1))) then
        ! Sending data to process i

        allocate(elements_send(size(enlist_b, 2)))
        nelements_send = 0
        do ele_B = 1, size(enlist_b, 2)
          ! Don't send non-owned elements
          if(ele_owner_b(ele_B) /= rank) cycle
          ! Don't send elements with non-intersecting bounding boxes
          if(.not. partition_bboxes_intersect(bbox(positions_b(:, enlist_b(:, ele_B))), parallel_bbox_a(:, :, i + 1))) cycle

          ! Mark the element for sending
          nelements_send = nelements_send + 1
          elements_send(nelements_send) = ele_B
        end do

        ! Build a set of nodes to send
        call allocate(nodes_send)
        do j = 1, nelements_send
          do k = 1, size(enlist_b, 1)
            call insert(nodes_send, enlist_b(k, elements_send(j)))
          end do
        end do        
        nnodes_send = key_count(nodes_send)

        allocate(send_enlist(nelements_send * size(enlist_b, 1)))

        ! Build an element-node graph for sending with contiguous indices
        call allocate(node_map)
        j = 0
        do k = 1, nnodes_send
          node = fetch(nodes_send, k)
          j = j + 1
          call insert(node_map, node, j)
        end do        
        do j = 1, nelements_send
          do k = 1, size(enlist_b, 1)
            send_enlist((j - 1) * size(enlist_b, 1) + k) = fetch(node_map, enlist_b(k, elements_send(j)))
          end do
        end do
        call deallocate(node_map)

        allocate(send_nodes_array(nnodes_send), send_positions(nnodes_send * size(positions_b, 1)))
        do j = 1, nnodes_send
          node = fetch(nodes_send, j)
          send_nodes_array(j) = node
          send_positions((j - 1) * size(positions_b, 1) + 1:j * size(positions_b, 1)) = positions_b(:, node)
        end do
        call deallocate(nodes_send)

      ! ### Mesh send buffer allocation ###
        if(nelements_send > 0) then
          call donor_ele_data(send_nodes_array, elements_send(:nelements_send), data)

          buffer_size = int_extent + int_extent +                          &
                    &   (size(send_enlist) * int_extent) + &
                    &   (size(send_positions)      * dp_extent ) + &
                    &   int_extent + size(data)

          allocate(send_buffer(i + 1)%v(buffer_size))
          position = 0

          call MPI_Pack(nelements_send, &
               & 1,                                                 &
               & MPI_INTEGER, send_buffer(i + 1)%v, buffer_size, position, mpi_comm, ierr)
          if(ierr /= MPI_SUCCESS) then
            libsupermesh_abort("MPI_Pack number of elements error")
          end if
          call MPI_Pack(nnodes_send, &
               & 1,                                                 &
               & MPI_INTEGER, send_buffer(i + 1)%v, buffer_size, position, mpi_comm, ierr)
          if(ierr /= MPI_SUCCESS) then
            libsupermesh_abort("MPI_Pack number of nodes error")
          end if
          call MPI_Pack(send_enlist,                &
               & size(send_enlist),                 &
               & MPI_INTEGER, send_buffer(i + 1)%v, buffer_size, position, mpi_comm, ierr)
          if(ierr /= MPI_SUCCESS) then
            libsupermesh_abort("MPI_Pack connectivity error")
          end if
          call MPI_Pack(send_positions,                     &
               & size(send_positions),                      &
               & MPI_DOUBLE_PRECISION, send_buffer(i + 1)%v, buffer_size, position, mpi_comm, ierr)
          if(ierr /= MPI_SUCCESS) then
            libsupermesh_abort("MPI_Pack values error")
          end if
          call MPI_Pack(size(data),                                 &
               & 1,                                                 &
               & MPI_INTEGER, send_buffer(i + 1)%v, buffer_size, position, mpi_comm, ierr)
          if(ierr /= MPI_SUCCESS) then
            libsupermesh_abort("MPI_Pack size of data error")
          end if
          call MPI_Pack(data,                                       &
               & size(data),                                        &
               & MPI_BYTE, send_buffer(i + 1)%v, buffer_size, position, mpi_comm, ierr)
          if(ierr /= MPI_SUCCESS) then
            libsupermesh_abort("MPI_Pack data error")
          end if

          deallocate(data)
        else
          buffer_size = int_extent + int_extent

          allocate(send_buffer(i + 1)%v(buffer_size))
          position = 0
          call MPI_Pack(nelements_send, &
               & 1,                                                 &
               & MPI_INTEGER, send_buffer(i + 1)%v, buffer_size, position, mpi_comm, ierr)
          if(ierr /= MPI_SUCCESS) then
            libsupermesh_abort("MPI_Pack number of elements error")
          end if
          call MPI_Pack(nnodes_send, &
               & 1,                                                 &
               & MPI_INTEGER, send_buffer(i + 1)%v, buffer_size, position, mpi_comm, ierr)
          if(ierr /= MPI_SUCCESS) then
            libsupermesh_abort("MPI_Pack number of nodes error")
          end if

        end if
        deallocate(send_nodes_array, send_enlist, send_positions)

#ifdef PROFILE
        if ( t0 .eq. -1 ) t0 = mpi_wtime()
#endif

        nsends = nsends + 1
        call MPI_Isend(send_buffer(i + 1)%v, buffer_size, MPI_PACKED, &
          & i, 0, mpi_comm, send_requests(nsends), ierr);  assert(ierr == MPI_SUCCESS)
    end if
    deallocate(elements_send)
  end do

  end subroutine step_2

  subroutine step_3(positions_a, enlist_a, ele_owner_a, &
                &   positions_b, enlist_b,       ele_owner_b, &
                &   unpack_data_b, intersection_calculation)
    ! dim x nnodes_a
    real, dimension(:, :), intent(in)       :: positions_a
    ! loc_a x nelements_a
    integer, dimension(:, :), intent(in)    :: enlist_a
    ! elements_a
    integer, dimension(:), intent(in)       :: ele_owner_a
    ! dim x nnodes_b
    real, dimension(:, :), intent(in)       :: positions_b
    ! loc_a x nelements_b
    integer, dimension(:, :), intent(in)    :: enlist_b
    ! elements_b
    integer, dimension(:), intent(in)       :: ele_owner_b

    integer                                 :: i, j, k, l, n_C, ele_A, ele_B, nintersections, ierr, position
    integer                                 :: nelements, nnodes
    real, dimension(:, :), allocatable      :: nodes_A, nodes_B
    real, dimension(:, :, :), allocatable   :: positions_c

    interface
      subroutine unpack_data_b(data_b)
        use iso_c_binding, only : c_int8_t
        implicit none
        integer(kind = c_int8_t), dimension(:), intent(in) :: data_b
      end subroutine unpack_data_b

      subroutine intersection_calculation(positions_a, positions_b, positions_c, nodes_b, ele_a, ele_b, local)
        use iso_c_binding, only : c_int8_t
        implicit none
        ! dim x loc_a
        real, dimension(:, :), intent(in) :: positions_a
        ! dim x loc_b
        real, dimension(:, :), intent(in) :: positions_b
        ! dim x loc_c x nelements_c
        real, dimension(:, :, :), intent(in) :: positions_c
        ! loc_b
        integer, dimension(:), intent(in) :: nodes_b
        integer, intent(in) :: ele_a
        integer, intent(in) :: ele_b
        logical, intent(in) :: local
      end subroutine intersection_calculation
    end interface

    integer :: buffer_size
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer, dimension(:), allocatable :: recv_enlist, statuses
    integer(kind = c_int8_t), dimension(:), allocatable :: data
    real, dimension(:), allocatable :: recv_positions
#ifdef OVERLAP_COMPUTE_COMMS
    integer(kind = c_int8_t), dimension(:), allocatable :: recv_buffer
#else
    integer(kind = c_int8_t), dimension(:), pointer :: recv_buffer
    type(buffer), dimension(:), allocatable :: recv_buffers
#endif

#ifndef OVERLAP_COMPUTE_COMMS
    allocate(recv_buffers(nprocs))
    do i = 0, nprocs - 1
      if(i == rank) cycle

      if(recv(i + 1)) then
        call MPI_Probe(i, MPI_ANY_TAG, mpi_comm, status, ierr);  assert(ierr == MPI_SUCCESS)
        call MPI_Get_count(status, MPI_PACKED, buffer_size, ierr);  assert(ierr == MPI_SUCCESS)

        allocate(recv_buffers(i + 1)%v(buffer_size))
        call MPI_Recv(recv_buffers(i + 1)%v, buffer_size, MPI_PACKED, & 
          & i, MPI_ANY_TAG, mpi_comm, status, ierr);  assert(ierr == MPI_SUCCESS)
      end if
    end do

    allocate(statuses(nsends * MPI_STATUS_SIZE))
    call MPI_Waitall(nsends, send_requests, statuses, ierr);  assert(ierr == MPI_SUCCESS)
    deallocate(statuses)
#ifdef PROFILE
    if ( t0 .gt. 0 ) then
      point_to_point = mpi_wtime() - t0
    else
      point_to_point = 0.0D0
    end if
#endif
#endif

    allocate(nodes_A(size(positions_a, 1), size(enlist_a, 1)), &
           & nodes_B(size(positions_a, 1), size(enlist_b, 1)), &
           & positions_c(size(positions_a, 1), size(positions_a, 1) + 1, &
                       & intersection_buffer_size(size(positions_a, 1), size(enlist_a, 1), size(enlist_b, 1))))
    nodes_A = 0.0
    nodes_B = 0.0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Parallel self-self runtime test !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PROFILE
    t1 = mpi_wtime()
#endif
    call rtree_intersection_finder_set_input(positions_a, enlist_a)
    do ele_B = 1, size(enlist_b, 2)
      if(ele_owner_b(ele_B) /= rank) cycle
      nodes_B = positions_b(:, enlist_b(:, ele_b))
      call rtree_intersection_finder_find(nodes_B)
      call rtree_intersection_finder_query_output(nintersections)
      do i = 1, nintersections
        call rtree_intersection_finder_get_output(ele_a, i)
        if(ele_owner_a(ele_A) /= rank) cycle
        nodes_A = positions_a(:, enlist_a(:, ele_A))

        call intersect_elements(nodes_A, nodes_B, positions_c, n_C)
        call intersection_calculation(nodes_A, nodes_B, positions_c(:, :, :n_C), enlist_b(:, ele_B), ele_A, ele_B, .true.)

      end do
    end do
#ifdef PROFILE
    local_compute = mpi_wtime() - t1
#endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Parallel self-other runtime test !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 0, nprocs - 1
      if(i == rank) cycle

      if(recv(i + 1)) then
#ifdef OVERLAP_COMPUTE_COMMS
        call MPI_Probe(i, MPI_ANY_TAG, mpi_comm, status, ierr);  assert(ierr == MPI_SUCCESS)
        call MPI_Get_count(status, MPI_PACKED, buffer_size, ierr);  assert(ierr == MPI_SUCCESS)

        allocate(recv_buffer(buffer_size))
        call MPI_Recv(recv_buffer, buffer_size, MPI_PACKED, & 
          & i, MPI_ANY_TAG, mpi_comm, status, ierr);  assert(ierr == MPI_SUCCESS)
#else
        recv_buffer => recv_buffers(i + 1)%v
        buffer_size = size(recv_buffer)
#endif
        position = 0
        call MPI_UnPack ( recv_buffer, buffer_size,  &
             & position, nelements,                        &
             & 1,                              &
             & MPI_INTEGER, mpi_comm, IERR)
        if(ierr /= MPI_SUCCESS) then
          libsupermesh_abort("MPI_UnPack number of elements error")
        end if
        allocate(recv_enlist(nelements * size(enlist_b, 1)))

        call MPI_UnPack ( recv_buffer, buffer_size,  &
             & position, nnodes,                    &
             & 1,                              &
             & MPI_INTEGER, mpi_comm, IERR)
        if(ierr /= MPI_SUCCESS) then
          libsupermesh_abort("MPI_UnPack number of nodes error")
        end if

        allocate(recv_positions(nnodes * (size(positions_b,1))))

        if(nelements <= 0) cycle

        call MPI_UnPack ( recv_buffer, buffer_size,  &
             & position, recv_enlist,     &
             & size(recv_enlist),         &
             & MPI_INTEGER, mpi_comm, IERR)
        if(ierr /= MPI_SUCCESS) then
          libsupermesh_abort("MPI_UnPack connectivity error")
        end if

        call MPI_UnPack ( recv_buffer, buffer_size,  &
             & position, recv_positions,          &
             & size(recv_positions),              &
             & MPI_DOUBLE_PRECISION, mpi_comm, IERR)
        if(ierr /= MPI_SUCCESS) then
          libsupermesh_abort("MPI_UnPack values error")
        end if

        call MPI_UnPack ( recv_buffer, buffer_size,  &
             & position, k,                               &
             & 1,                                         &
             & MPI_INTEGER, mpi_comm, IERR)
        if(ierr /= MPI_SUCCESS) then
          libsupermesh_abort("MPI_UnPack size of data error")
        end if

        allocate(data(k))
        call MPI_UnPack ( recv_buffer, buffer_size,  &
             & position, data,                           &
             & k,                                         &
             & MPI_BYTE, mpi_comm, IERR)
        if(ierr /= MPI_SUCCESS) then
          libsupermesh_abort("MPI_UnPack data error")
        end if
#ifdef OVERLAP_COMPUTE_COMMS
        deallocate(recv_buffer)
#endif

        call unpack_data_b(data)

#ifdef PROFILE
    t1 = mpi_wtime()
#endif
        do ele_B = 1, nelements
          do k = 1, size(enlist_b, 1)
            do j = 1, size(positions_b, 1)
              l = recv_enlist((ele_B - 1) * size(enlist_b, 1) + k)            
              nodes_B(j, k) = recv_positions((l - 1) * size(positions_b, 1) + j)
            end do
          end do

          call rtree_intersection_finder_find(nodes_B)
          call rtree_intersection_finder_query_output(nintersections)
          do k = 1, nintersections
            call rtree_intersection_finder_get_output(ele_A, k)
            if(ele_owner_a(ele_A) /= rank) cycle
            nodes_A = positions_a(:, enlist_a(:, ele_A))

            call intersect_elements(nodes_A, nodes_B, positions_c, n_C)
            call intersection_calculation(nodes_A, nodes_B, positions_c(:, :, :n_C), recv_enlist((ele_B - 1) * size(enlist_b, 1) + 1:ele_B * size(enlist_B, 1)), ele_A, ele_B, .false.)
          end do
        end do
#ifdef PROFILE
    remote_compute = mpi_wtime() - t1
#endif
        deallocate(recv_enlist, recv_positions, data)
      end if
    end do

    deallocate(positions_c, nodes_A, nodes_B)

#ifdef OVERLAP_COMPUTE_COMMS
    allocate(statuses(nsends * MPI_STATUS_SIZE))
    call MPI_Waitall(nsends, send_requests, statuses, ierr);  assert(ierr == MPI_SUCCESS)
    deallocate(statuses)
#else
    do i = 1, nprocs
      if(associated(recv_buffers(i)%v)) deallocate(recv_buffers(i)%v)
    end do
    deallocate(recv_buffers)
#endif

  end subroutine step_3

  subroutine initialise_parallel_supermesh(comm)
    integer, optional, intent(in) :: comm

    integer  :: ierr
    
    if(parallel_supermesh_allocated) call finalise_parallel_supermesh()
    
    if(present(comm)) then
      mpi_comm = comm
    else
      mpi_comm = MPI_COMM_WORLD
    end if

    call MPI_Comm_rank(mpi_comm, rank, ierr);  assert(ierr == MPI_SUCCESS)
    call MPI_Comm_size(mpi_comm, nprocs, ierr);  assert(ierr == MPI_SUCCESS)

    allocate(send_buffer(nprocs), &
           & send_requests(nprocs), &
           & recv(nprocs))
    send_requests = MPI_REQUEST_NULL
    
    parallel_supermesh_allocated = .true.
    
  end subroutine initialise_parallel_supermesh

  subroutine finalise_parallel_supermesh()
    integer :: i

    if(.not. parallel_supermesh_allocated) return

    do i = 1, nprocs
      if(associated(send_buffer(i)%v)) deallocate(send_buffer(i)%v)
    end do
    deallocate(send_buffer, &
             & send_requests, &
             & recv)
    if(allocated(bbox_a)) deallocate(bbox_a, bbox_b, parallel_bbox_a, parallel_bbox_b)

    parallel_supermesh_allocated = .false.

  end subroutine finalise_parallel_supermesh

  subroutine print_times(comm, unit)
    integer, optional, intent(in) :: comm
    integer, optional, intent(in) :: unit
    
#ifdef PROFILE
    integer :: lcomm, lunit
  
    integer :: ierr, nprocs
    real :: point_to_point_max, point_to_point_min, point_to_point_sum
    real :: all_to_all_max, all_to_all_min, all_to_all_sum
    real :: local_compute_max, local_compute_min, local_compute_sum
    real :: remote_compute_max, remote_compute_min, remote_compute_sum
    
    if(present(comm)) then
      lcomm = comm
    else
      lcomm = MPI_COMM_WORLD
    end if
    if(present(unit)) then
      lunit = unit
    else
      lunit = debug_log_unit
    end if
    
    call MPI_Comm_size(lcomm, nprocs, ierr);  assert(ierr == MPI_SUCCESS)

    call MPI_Allreduce(point_to_point, point_to_point_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, lcomm, ierr);  assert(ierr == MPI_SUCCESS)
    call MPI_Allreduce(point_to_point, point_to_point_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, lcomm, ierr);  assert(ierr == MPI_SUCCESS)
    call MPI_Allreduce(point_to_point, point_to_point_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, lcomm, ierr);  assert(ierr == MPI_SUCCESS)
    call MPI_Allreduce(all_to_all,     all_to_all_min,     1, MPI_DOUBLE_PRECISION, MPI_MIN, lcomm, ierr);  assert(ierr == MPI_SUCCESS)
    call MPI_Allreduce(all_to_all,     all_to_all_max,     1, MPI_DOUBLE_PRECISION, MPI_MAX, lcomm, ierr);  assert(ierr == MPI_SUCCESS)
    call MPI_Allreduce(all_to_all,     all_to_all_sum,     1, MPI_DOUBLE_PRECISION, MPI_SUM, lcomm, ierr);  assert(ierr == MPI_SUCCESS)
    call MPI_Allreduce(local_compute,  local_compute_min,  1, MPI_DOUBLE_PRECISION, MPI_MIN, lcomm, ierr);  assert(ierr == MPI_SUCCESS)
    call MPI_Allreduce(local_compute,  local_compute_max,  1, MPI_DOUBLE_PRECISION, MPI_MAX, lcomm, ierr);  assert(ierr == MPI_SUCCESS)
    call MPI_Allreduce(local_compute,  local_compute_sum,  1, MPI_DOUBLE_PRECISION, MPI_SUM, lcomm, ierr);  assert(ierr == MPI_SUCCESS)
    call MPI_Allreduce(remote_compute, remote_compute_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, lcomm, ierr);  assert(ierr == MPI_SUCCESS)      
    call MPI_Allreduce(remote_compute, remote_compute_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, lcomm, ierr);  assert(ierr == MPI_SUCCESS)
    call MPI_Allreduce(remote_compute, remote_compute_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, lcomm, ierr);  assert(ierr == MPI_SUCCESS)

    if(rank == 0) then
      write(unit, "(a,e26.18e3,a,e26.18e3,a,e26.18e3)") "All to all time, min, max, average     = ", all_to_all_min,     ", ", all_to_all_max,     ", ", all_to_all_sum     / nprocs
      write(unit, "(a,e26.18e3,a,e26.18e3,a,e26.18e3)") "Point to point time, min, max, average = ", point_to_point_min, ", ", point_to_point_max, ", ", point_to_point_sum / nprocs
      write(unit, "(a,e26.18e3,a,e26.18e3,a,e26.18e3)") "Local compute time, min, max, average  = ", local_compute_min,  ", ", local_compute_max,  ", ", local_compute_sum  / nprocs
      write(unit, "(a,e26.18e3,a,e26.18e3,a,e26.18e3)") "Remote compute time, min, max, average = ", remote_compute_min, ", ", remote_compute_max, ", ", remote_compute_sum / nprocs
    end if
#endif

  end subroutine print_times

  pure function bbox(coords)
    ! dim x loc
    real, dimension(:, :), intent(in) :: coords

    real, dimension(2, size(coords, 1)) :: bbox

    integer :: i, j

    bbox(1, :) = coords(:, 1)
    bbox(2, :) = coords(:, 1)
    do i = 2, size(coords, 2)
      do j = 1, size(coords, 1)
        bbox(1, j) = min(bbox(1, j), coords(j, i))
        bbox(2, j) = max(bbox(2, j), coords(j, i))
      end do
    end do

  end function bbox

  pure function partition_bbox_real(positions, enlist, ele_owner, rank) result(bbox)
    ! dim x nnodes
    real, dimension(:, :), intent(in) :: positions
    ! loc x nelements
    integer, dimension(:, :), intent(in) :: enlist
    ! nelements
    integer, dimension(:), intent(in) :: ele_owner
    integer, intent(in) :: rank
    
    real, dimension(2, size(positions, 1)) :: bbox

    integer :: dim, ele, ele_0, i, nelements
    real, dimension(size(positions, 1), size(enlist, 1)) :: element
    
    dim = size(positions, 1)
    nelements = size(enlist, 2)    
    if(nelements == 0) then
      bbox = huge(0.0)
      return
    end if
    
    ele_0 = 1
    do while(ele_owner(ele_0) /= rank)
      ele_0 = ele_0 + 1
      if(ele_0 > nelements) then
        bbox = huge(0.0)
        return
      end if
    end do    
    element = positions(:, enlist(:, ele_0))
    do i = 1, dim
      bbox(1, i) = minval(element(1, :))
      bbox(2, i) = maxval(element(2, :))
    end do
    
    do ele = ele_0 + 1, nelements
      if(ele_owner(ele) /= rank) cycle
      element = positions(:, enlist(:, ele))
      do i = 1, dim
        bbox(1, i) = min(bbox(1, i), minval(element(i, :)))
        bbox(2, i) = max(bbox(2, i), maxval(element(i, :)))
      end do
    end do

  end function partition_bbox_real
  
  pure function partition_bbox_vector_field(positions, ele_owner, rank) result(bbox)
    use libsupermesh_fields, only : vector_field
    type(vector_field), intent(in) :: positions
    ! nelements
    integer, dimension(:), intent(in) :: ele_owner
    integer, intent(in) :: rank
    
    real, dimension(2, size(positions%val, 1)) :: bbox
    
    bbox = partition_bbox(positions%val, reshape(positions%mesh%ndglno, (/positions%mesh%loc, positions%mesh%nelements/)), ele_owner, rank)
  
  end function partition_bbox_vector_field

  pure function partition_bboxes_intersect(bbox_1, bbox_2) result(intersect)
    ! 2 x dim
    real, dimension(:, :), intent(in) :: bbox_1
    ! 2 x dim
    real, dimension(:, :), intent(in) :: bbox_2

    logical :: intersect

    integer :: i
    
    do i = 1, size(bbox_1, 2)
      if(bbox_2(2, i) <= bbox_1(1, i) .or. bbox_2(1, i) >= bbox_1(2, i)) then
        intersect = .false.
        return
      end if
    end do
    intersect = .true.

  end function partition_bboxes_intersect

end module libsupermesh_parallel_supermesh
