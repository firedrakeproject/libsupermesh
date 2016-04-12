!  Copyright (C) 2016 The University of Edinburgh
!
!  The file is part of libsupermesh
!    
!  This library is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation;
!  version 2.1 of the License.
!
!  This library is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public
!  License along with this library; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

#include "libsupermesh_debug.h"

module libsupermesh_parallel_supermesh

  use iso_c_binding, only : c_int8_t
  
  use libsupermesh_debug, only : abort_pinpoint
  use libsupermesh_debug_parameters, only : debug_log_unit
  use libsupermesh_integer_hash_table, only : integer_hash_table, allocate, &
    & insert, fetch, deallocate
  use libsupermesh_integer_set, only : integer_set, allocate, deallocate, &
    & insert, key_count, fetch
  use libsupermesh_intersection_finder, only : rtree_type, allocate, query, &
    & deallocate
  use libsupermesh_supermesh, only : max_n_simplices_c, intersect_elements

  implicit none

#include <mpif.h>

  private

  public :: parallel_supermesh, print_times
  public :: partition_bbox, partition_bboxes_intersect

  type buffer
    integer(kind = c_int8_t), dimension(:), pointer :: v => null()
  end type buffer

  logical, save :: parallel_supermesh_allocated = .false.
  
  real, dimension(:, :), allocatable, save :: bbox_a, bbox_b
  real, dimension(:, :, :), allocatable, save :: parallel_bbox_a, &
    & parallel_bbox_b

  integer :: nsends
  integer, dimension(:), allocatable, save :: send_requests
  logical, dimension(:), allocatable, save :: recv
  type(buffer), dimension(:), allocatable, save :: send_buffer

  integer, save :: mpi_comm, nprocs, rank
#ifdef PROFILE
  real :: t_0, all_to_all_time, local_compute_time, remote_compute_time
#ifndef OVERLAP_COMPUTE_COMMS
  real :: point_to_point_time
#endif
#endif

contains

  subroutine parallel_supermesh(positions_a, enlist_a, ele_owner_a, &
                             &  positions_b, enlist_b, ele_owner_b, &
                             &  pack_data_b, unpack_data_b, intersection_calculation, &
                             &  comm)
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
    integer, optional, intent(in) :: comm
    interface
      subroutine pack_data_b(nodes_b, eles_b, data_b)
        use iso_c_binding, only : c_int8_t
        implicit none
        integer, dimension(:), intent(in) :: nodes_b
        integer, dimension(:), intent(in) :: eles_b
        integer(kind = c_int8_t), dimension(:), allocatable, intent(out) :: data_b
      end subroutine pack_data_b

      subroutine unpack_data_b(nnodes_b, nelements_b, data_b)
        use iso_c_binding, only : c_int8_t
        implicit none
        integer, intent(in) :: nnodes_b
        integer, intent(in) :: nelements_b
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

    ! 0. Allocate and initialise arrays
    call initialise_parallel_supermesh(comm)

    ! 1. Calculate and communicate bounding box data for donor and target
    call step_1(positions_a, enlist_a, ele_owner_a, positions_b, enlist_b, &
      & ele_owner_b)

    ! 2. Use bounding box data to cull donor mesh
    call step_2(positions_b, enlist_b, ele_owner_b, pack_data_b)

    ! 5. Supermesh and call user specified element unpack and calculation functions
    call step_3(positions_a, enlist_a, ele_owner_a, positions_b, enlist_b, &
      & ele_owner_b, unpack_data_b, intersection_calculation)

    ! 6. Deallocate
    call finalise_parallel_supermesh()

  end subroutine parallel_supermesh

  subroutine step_1(positions_a, enlist_a, ele_owner_a, &
                  & positions_b, enlist_b, ele_owner_b)
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
#ifdef PROFILE
    real :: t_1
#endif

    dim = size(positions_a, 1)
    
    allocate(bbox_a(2, dim))
    bbox_a = partition_bbox(positions_a, enlist_a, ele_owner_a, rank)

    allocate(bbox_b(2, dim))
    bbox_b = partition_bbox(positions_b, enlist_b, ele_owner_b, rank)

    allocate(parallel_bbox_a(2, dim, nprocs))
    allocate(parallel_bbox_b(2, dim, nprocs))

#ifdef PROFILE
    t_1 = mpi_wtime()
#endif
    call MPI_Allgather(bbox_a, 2 * dim, MPI_DOUBLE_PRECISION, &
      & parallel_bbox_a, 2 * dim, MPI_DOUBLE_PRECISION, &
      & mpi_comm, ierr);  assert(ierr == MPI_SUCCESS)

    call MPI_Allgather(bbox_b, 2 * dim, MPI_DOUBLE_PRECISION, &
      & parallel_bbox_b, 2 * dim, MPI_DOUBLE_PRECISION, &
      & mpi_comm, ierr);  assert(ierr == MPI_SUCCESS)
#ifdef PROFILE
    all_to_all_time = mpi_wtime() - t_1
#endif

  end subroutine step_1

  subroutine step_2(positions_b, enlist_b, ele_owner_b, pack_data_b)
    ! dim x nnodes_b
    real, dimension(:, :), intent(in) :: positions_b
    ! loc_b x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    ! elements_b
    integer, dimension(:), intent(in) :: ele_owner_b

    interface
      subroutine pack_data_b(nodes_b, eles_b, data_b)
        use iso_c_binding, only : c_int8_t
        implicit none
        integer, dimension(:), intent(in) :: nodes_b
        integer, dimension(:), intent(in) :: eles_b
        integer(kind = c_int8_t), dimension(:), allocatable, intent(out) :: data_b
      end subroutine pack_data_b
    end interface

    integer :: buffer_size, dim, dp_extent, ele_b, i, ierr, int_extent, j, k, &
      & loc_b, nelements_b, nelements_send, nnodes_send, node, position
    integer, dimension(:), allocatable :: elements_send, send_enlist, &
      & send_nodes_array
    integer(kind = c_int8_t), dimension(:), allocatable :: data
    real, dimension(:), allocatable :: send_positions
    type(integer_hash_table) :: node_map
    type(integer_set) :: nodes_send

#ifdef PROFILE
    t_0 = -1.0D0
#endif

    dim = size(positions_b, 1)
    loc_b = size(enlist_b, 1)
    nelements_b = size(enlist_b, 2)

    call MPI_Type_extent(MPI_DOUBLE_PRECISION, dp_extent, ierr);  assert(ierr == MPI_SUCCESS)
    call MPI_Type_extent(MPI_INTEGER, int_extent, ierr);  assert(ierr == MPI_SUCCESS)

    allocate(elements_send(nelements_b))
    nsends = 0
    do i = 0, nprocs - 1
      if(i == rank) then
        recv(i + 1) = .false.
        cycle
      end if

      recv(i + 1) = partition_bboxes_intersect(bbox_a, parallel_bbox_b(:, :, i + 1))

      if(partition_bboxes_intersect(bbox_b, parallel_bbox_a(:, :, i + 1))) then
        nelements_send = 0
        do ele_b = 1, nelements_b
          ! Don't send non-owned elements
          if(ele_owner_b(ele_b) /= rank) cycle
          ! Don't send elements with non-intersecting bounding boxes
          if(.not. partition_bboxes_intersect(bbox(positions_b(:, enlist_b(:, ele_b))), parallel_bbox_a(:, :, i + 1))) cycle

          ! Mark the element for sending
          nelements_send = nelements_send + 1
          elements_send(nelements_send) = ele_b
        end do

        ! Build a set of nodes to send
        call allocate(nodes_send)
        do j = 1, nelements_send
          do k = 1, loc_b
            call insert(nodes_send, enlist_b(k, elements_send(j)))
          end do
        end do        
        nnodes_send = key_count(nodes_send)

        ! Build an element-node graph for sending with contiguous indices
        call allocate(node_map)
        j = 0
        do k = 1, nnodes_send
          node = fetch(nodes_send, k)
          j = j + 1
          call insert(node_map, node, j)
        end do        
        allocate(send_enlist(loc_b * nelements_send))
        do j = 1, nelements_send
          do k = 1, loc_b
            send_enlist((j - 1) * loc_b + k) = fetch(node_map, enlist_b(k, elements_send(j)))
          end do
        end do
        call deallocate(node_map)

        allocate(send_nodes_array(nnodes_send), send_positions(dim * nnodes_send))
        do j = 1, nnodes_send
          node = fetch(nodes_send, j)
          send_nodes_array(j) = node
          send_positions((j - 1) * dim + 1:j * dim) = positions_b(:, node)
        end do
        call deallocate(nodes_send)

        if(nelements_send > 0) then
          call pack_data_b(send_nodes_array, elements_send(:nelements_send), data)

          buffer_size = int_extent + int_extent + &
                      & (size(send_enlist) * int_extent) + &
                      & (size(send_positions) * dp_extent) + &
                      & int_extent + size(data)
          allocate(send_buffer(i + 1)%v(buffer_size))
          position = 0
          call MPI_Pack(nelements_send, 1, MPI_INTEGER, &
            & send_buffer(i + 1)%v, buffer_size, position, mpi_comm, ierr);  assert(ierr == MPI_SUCCESS)
          call MPI_Pack(nnodes_send, 1, MPI_INTEGER, &
            & send_buffer(i + 1)%v, buffer_size, position, mpi_comm, ierr);  assert(ierr == MPI_SUCCESS)
          call MPI_Pack(send_enlist, size(send_enlist), MPI_INTEGER, &
            & send_buffer(i + 1)%v, buffer_size, position, mpi_comm, ierr);  assert(ierr == MPI_SUCCESS)
          call MPI_Pack(send_positions, size(send_positions), MPI_DOUBLE_PRECISION, &
            & send_buffer(i + 1)%v, buffer_size, position, mpi_comm, ierr);  assert(ierr == MPI_SUCCESS)
          call MPI_Pack(size(data), 1, MPI_INTEGER, &
            & send_buffer(i + 1)%v, buffer_size, position, mpi_comm, ierr);  assert(ierr == MPI_SUCCESS)
          call MPI_Pack(data, size(data), MPI_BYTE, &
            & send_buffer(i + 1)%v, buffer_size, position, mpi_comm, ierr);  assert(ierr == MPI_SUCCESS)

          deallocate(data)
        else
          buffer_size = int_extent + int_extent
          allocate(send_buffer(i + 1)%v(buffer_size))
          position = 0
          call MPI_Pack(nelements_send, 1, MPI_INTEGER, &
            & send_buffer(i + 1)%v, buffer_size, position, mpi_comm, ierr);  assert(ierr == MPI_SUCCESS)
          call MPI_Pack(nnodes_send, 1, MPI_INTEGER, &
            & send_buffer(i + 1)%v, buffer_size, position, mpi_comm, ierr);  assert(ierr == MPI_SUCCESS)
        end if
        deallocate(send_enlist, send_nodes_array, send_positions)

#ifdef PROFILE
        if(t_0 < 0.0D0) t_0 = mpi_wtime()
#endif

        nsends = nsends + 1
        call MPI_Isend(send_buffer(i + 1)%v, buffer_size, MPI_PACKED, &
          & i, 0, mpi_comm, send_requests(nsends), ierr);  assert(ierr == MPI_SUCCESS)
    end if
  end do
  deallocate(elements_send)

  end subroutine step_2

  subroutine step_3(positions_a, enlist_a, ele_owner_a, &
                  & positions_b, enlist_b, ele_owner_b, &
                  & unpack_data_b, intersection_calculation)
    ! dim x nnodes_a
    real, dimension(:, :), intent(in) :: positions_a
    ! loc_a x nelements_a
    integer, dimension(:, :), intent(in) :: enlist_a
    ! elements_a
    integer, dimension(:), intent(in) :: ele_owner_a
    ! dim x nnodes_b
    real, dimension(:, :), intent(in) :: positions_b
    ! loc_a x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    ! elements_b
    integer, dimension(:), intent(in) :: ele_owner_b

    interface
      subroutine unpack_data_b(nnodes_b, nelements_b, data_b)
        use iso_c_binding, only : c_int8_t
        implicit none
        integer, intent(in) :: nnodes_b
        integer, intent(in) :: nelements_b
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

    integer :: buffer_size, dim, ele_a, ele_b, i, ierr, j, k, l, loc_a, loc_b, &
      & n_c, ndata, nelements, nelements_b, nnodes, position
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer, dimension(:), allocatable :: eles_a, recv_enlist, statuses
    integer(kind = c_int8_t), dimension(:), allocatable :: data
    real, dimension(:), allocatable :: recv_positions
    real, dimension(:, :), allocatable :: nodes_a, nodes_b
    real, dimension(:, :, :), allocatable :: positions_c
    type(rtree_type) :: rtree
#ifdef OVERLAP_COMPUTE_COMMS
    integer(kind = c_int8_t), dimension(:), allocatable :: recv_buffer
#else
    integer(kind = c_int8_t), dimension(:), pointer :: recv_buffer
    type(buffer), dimension(:), allocatable :: recv_buffers
#endif
#ifdef PROFILE
    real :: t_1
#endif

#ifndef OVERLAP_COMPUTE_COMMS
    allocate(recv_buffers(nprocs))
    do i = 0, nprocs - 1
      if(.not. recv(i + 1)) cycle
      
      call MPI_Probe(i, MPI_ANY_TAG, mpi_comm, status, ierr);  assert(ierr == MPI_SUCCESS)
      call MPI_Get_count(status, MPI_PACKED, buffer_size, ierr);  assert(ierr == MPI_SUCCESS)

      allocate(recv_buffers(i + 1)%v(buffer_size))
      call MPI_Recv(recv_buffers(i + 1)%v, buffer_size, MPI_PACKED, & 
        & i, MPI_ANY_TAG, mpi_comm, status, ierr);  assert(ierr == MPI_SUCCESS)
    end do

    allocate(statuses(nsends * MPI_STATUS_SIZE))
    call MPI_Waitall(nsends, send_requests, statuses, ierr);  assert(ierr == MPI_SUCCESS)
    deallocate(statuses)
    do i = 1, nprocs
      if(associated(send_buffer(i)%v)) then
        deallocate(send_buffer(i)%v)
        nullify(send_buffer(i)%v)
      end if
    end do
#ifdef PROFILE
    if(t_0 >= 0.0D0) then
      point_to_point_time = mpi_wtime() - t_0
    else
      point_to_point_time = 0.0D0
    end if
#endif
#endif

    dim = size(positions_a, 1)
    loc_a = size(enlist_a, 1)
    loc_b = size(enlist_b, 1)
    nelements_b = size(enlist_b, 2)

    allocate(nodes_a(dim, loc_a), nodes_b(dim, loc_b), &
           & positions_c(dim, dim + 1, max_n_simplices_c(dim, loc_a, loc_b)))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Self-self integration !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PROFILE
    t_1 = mpi_wtime()
#endif
    call allocate(rtree, positions_a, enlist_a)
    do ele_b = 1, nelements_b
      if(ele_owner_b(ele_b) /= rank) cycle
      
      nodes_b = positions_b(:, enlist_b(:, ele_b))
      call query(rtree, nodes_b, eles_a)
      do i = 1, size(eles_a)
        ele_a = eles_a(i)
        if(ele_owner_a(ele_a) /= rank) cycle
        
        nodes_a = positions_a(:, enlist_a(:, ele_a))
        call intersect_elements(nodes_a, nodes_b, positions_c, n_c)
        call intersection_calculation(nodes_a, nodes_b, &
          & positions_c(:, :, :n_c), enlist_b(:, ele_b), ele_a, ele_b, &
          & local = .true.)
      end do
      deallocate(eles_a)
    end do
#ifdef PROFILE
    local_compute_time = mpi_wtime() - t_1
#endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Self-received integration !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 0, nprocs - 1
      if(.not. recv(i + 1)) cycle
    
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
      call MPI_Unpack(recv_buffer, buffer_size, position, nelements, &
        & 1, MPI_INTEGER, mpi_comm, ierr);  assert(ierr == MPI_SUCCESS)

      if(nelements <= 0) then
#ifdef OVERLAP_COMPUTE_COMMS
        deallocate(recv_buffer)
#else
        deallocate(recv_buffers(i + 1)%v)
        nullify(recv_buffers(i + 1)%v)
#endif
        cycle
      end if

      call MPI_Unpack(recv_buffer, buffer_size, position, nnodes, &
        & 1, MPI_INTEGER, mpi_comm, ierr);  assert(ierr == MPI_SUCCESS)
      allocate(recv_enlist(loc_b * nelements), &
             & recv_positions(dim * nnodes))
      call MPI_Unpack(recv_buffer, buffer_size, position, recv_enlist, &
        & size(recv_enlist), MPI_INTEGER, mpi_comm, ierr);  assert(ierr == MPI_SUCCESS)
      call MPI_Unpack(recv_buffer, buffer_size, position, recv_positions, &
        & size(recv_positions), MPI_DOUBLE_PRECISION, mpi_comm, ierr);  assert(ierr == MPI_SUCCESS)
      call MPI_Unpack(recv_buffer, buffer_size, position, ndata, &
        & 1, MPI_INTEGER, mpi_comm, ierr);  assert(ierr == MPI_SUCCESS)
      allocate(data(ndata))
      call MPI_Unpack(recv_buffer, buffer_size, position, data, &
        & ndata, MPI_BYTE, mpi_comm, ierr);  assert(ierr == MPI_SUCCESS)
        
#ifdef OVERLAP_COMPUTE_COMMS
      deallocate(recv_buffer)
#else
      deallocate(recv_buffers(i + 1)%v)
      nullify(recv_buffers(i + 1)%v)
#endif

      call unpack_data_b(nnodes, nelements, data)

#ifdef PROFILE
      t_1 = mpi_wtime()
#endif
      do ele_b = 1, nelements
        do k = 1, loc_b
          do j = 1, dim
            l = recv_enlist((ele_b - 1) * loc_b + k)            
            nodes_b(j, k) = recv_positions((l - 1) * dim + j)
          end do
        end do

        call query(rtree, nodes_b, eles_a)
        do k = 1, size(eles_a)
          ele_a = eles_a(k)
          if(ele_owner_a(ele_a) /= rank) cycle
          
          nodes_a = positions_a(:, enlist_a(:, ele_a))
          call intersect_elements(nodes_a, nodes_b, positions_c, n_c)
          call intersection_calculation(nodes_a, nodes_b, &
            & positions_c(:, :, :n_c), &
            & recv_enlist((ele_b - 1) * loc_b + 1:ele_b * loc_b), ele_a, &
            & ele_b, local = .false.)
        end do
        deallocate(eles_a)
      end do
#ifdef PROFILE
      remote_compute_time = mpi_wtime() - t_1
#endif

      deallocate(recv_enlist, recv_positions, data)
    end do

    call deallocate(rtree)
    deallocate(nodes_a, nodes_b, positions_c)

#ifdef OVERLAP_COMPUTE_COMMS
    allocate(statuses(nsends * MPI_STATUS_SIZE))
    call MPI_Waitall(nsends, send_requests, statuses, ierr);  assert(ierr == MPI_SUCCESS)
    deallocate(statuses)
    do i = 1, nprocs
      if(associated(send_buffer(i)%v)) then
        deallocate(send_buffer(i)%v)
        nullify(send_buffer(i)%v)
      end if
    end do
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

    allocate(send_buffer(nprocs), send_requests(nprocs), recv(nprocs))
    send_requests = MPI_REQUEST_NULL
    
#ifdef PROFILE
    t_0 = -1.0D0
    all_to_all_time = 0.0D0
#ifndef OVERLAP_COMPUTE_COMMS
    point_to_point_time = 0.0D0
#endif
    local_compute_time = 0.0D0
    remote_compute_time = 0.0D0
#endif
    
    parallel_supermesh_allocated = .true.
    
  end subroutine initialise_parallel_supermesh

  subroutine finalise_parallel_supermesh()
    integer :: i

    if(.not. parallel_supermesh_allocated) return

    do i = 1, nprocs
      if(associated(send_buffer(i)%v)) deallocate(send_buffer(i)%v)
    end do
    deallocate(send_buffer, send_requests, recv)
    if(allocated(bbox_a)) deallocate(bbox_a, bbox_b, parallel_bbox_a, parallel_bbox_b)

    parallel_supermesh_allocated = .false.

  end subroutine finalise_parallel_supermesh

  subroutine print_times(comm, unit)
    integer, optional, intent(in) :: comm
    integer, optional, intent(in) :: unit
    
#ifdef PROFILE
    integer :: lcomm, lunit
  
    integer :: ierr, nprocs, rank
    real :: all_to_all_max, all_to_all_min, all_to_all_sum
#ifndef OVERLAP_COMPUTE_COMMS
    real :: point_to_point_max, point_to_point_min, point_to_point_sum
#endif
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
    
    call MPI_Comm_rank(lcomm, rank, ierr);  assert(ierr == MPI_SUCCESS)
    call MPI_Comm_size(lcomm, nprocs, ierr);  assert(ierr == MPI_SUCCESS)

    call MPI_Allreduce(all_to_all_time,     all_to_all_min,     1, MPI_DOUBLE_PRECISION, MPI_MIN, lcomm, ierr);  assert(ierr == MPI_SUCCESS)
    call MPI_Allreduce(all_to_all_time,     all_to_all_max,     1, MPI_DOUBLE_PRECISION, MPI_MAX, lcomm, ierr);  assert(ierr == MPI_SUCCESS)
    call MPI_Allreduce(all_to_all_time,     all_to_all_sum,     1, MPI_DOUBLE_PRECISION, MPI_SUM, lcomm, ierr);  assert(ierr == MPI_SUCCESS)
#ifndef OVERLAP_COMPUTE_COMMS
    call MPI_Allreduce(point_to_point_time, point_to_point_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, lcomm, ierr);  assert(ierr == MPI_SUCCESS)
    call MPI_Allreduce(point_to_point_time, point_to_point_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, lcomm, ierr);  assert(ierr == MPI_SUCCESS)
    call MPI_Allreduce(point_to_point_time, point_to_point_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, lcomm, ierr);  assert(ierr == MPI_SUCCESS)
#endif
    call MPI_Allreduce(local_compute_time,  local_compute_min,  1, MPI_DOUBLE_PRECISION, MPI_MIN, lcomm, ierr);  assert(ierr == MPI_SUCCESS)
    call MPI_Allreduce(local_compute_time,  local_compute_max,  1, MPI_DOUBLE_PRECISION, MPI_MAX, lcomm, ierr);  assert(ierr == MPI_SUCCESS)
    call MPI_Allreduce(local_compute_time,  local_compute_sum,  1, MPI_DOUBLE_PRECISION, MPI_SUM, lcomm, ierr);  assert(ierr == MPI_SUCCESS)
    call MPI_Allreduce(remote_compute_time, remote_compute_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, lcomm, ierr);  assert(ierr == MPI_SUCCESS)      
    call MPI_Allreduce(remote_compute_time, remote_compute_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, lcomm, ierr);  assert(ierr == MPI_SUCCESS)
    call MPI_Allreduce(remote_compute_time, remote_compute_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, lcomm, ierr);  assert(ierr == MPI_SUCCESS)

    if(rank == 0) then
      write(lunit, "(a,e26.18e3,a,e26.18e3,a,e26.18e3)") "All to all time, min, max, mean     = ", all_to_all_min,     ", ", all_to_all_max,     ", ", all_to_all_sum     / nprocs
#ifndef OVERLAP_COMPUTE_COMMS
      write(lunit, "(a,e26.18e3,a,e26.18e3,a,e26.18e3)") "Point to point time, min, max, mean = ", point_to_point_min, ", ", point_to_point_max, ", ", point_to_point_sum / nprocs
#endif
      write(lunit, "(a,e26.18e3,a,e26.18e3,a,e26.18e3)") "Local compute time, min, max, mean  = ", local_compute_min,  ", ", local_compute_max,  ", ", local_compute_sum  / nprocs
      write(lunit, "(a,e26.18e3,a,e26.18e3,a,e26.18e3)") "Remote compute time, min, max, mean = ", remote_compute_min, ", ", remote_compute_max, ", ", remote_compute_sum / nprocs
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

  pure function partition_bbox(positions, enlist, ele_owner, rank) result(bbox)
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
      bbox = huge(0.0D0)
      return
    end if
    
    ele_0 = 1
    do while(ele_owner(ele_0) /= rank)
      ele_0 = ele_0 + 1
      if(ele_0 > nelements) then
        bbox = huge(0.0D0)
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

  end function partition_bbox

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
