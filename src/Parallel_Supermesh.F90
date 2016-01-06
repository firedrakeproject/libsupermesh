#include "fdebug.h"

module libsupermesh_parallel_supermesh

  use iso_fortran_env, only : output_unit
  use iso_c_binding, only : c_int8_t
  use libsupermesh_integer_hash_table
  use libsupermesh_fldebug
  use libsupermesh_read_halos
  use libsupermesh_intersection_finder
  use libsupermesh_tri_intersection_module
  use libsupermesh_tet_intersection_module
  use libsupermesh_construction
  use libsupermesh_fields, only : triangle_area
  use libsupermesh_integer_set

  implicit none

  private

  public :: parallel_supermesh

  type pointer_real
    real, dimension(:), pointer :: p
  end type pointer_real

  type pointer_byte
    integer(kind = c_int8_t), dimension(:), pointer :: p
  end type pointer_byte

  type pointer_integer
    integer, dimension(:), pointer :: p
  end type pointer_integer

  logical, save :: parallel_supermesh_allocated = .false.
  real, dimension(:,:), allocatable, save :: bbox_a, bbox_b
  real, dimension(:,:,:), allocatable, save :: parallel_bbox_a, parallel_bbox_b

  type(pointer_integer), dimension(:), allocatable, save :: send_element_uns

  type(pointer_byte), dimension(:), allocatable, save :: recv_buffer, send_buffer
  integer, dimension(:,:), allocatable, save :: number_of_elements_and_nodes_to_receive, number_of_elements_and_nodes_to_send
  type(pointer_integer), dimension(:), allocatable, save :: send_nodes_connectivity, recv_nodes_connectivity
  type(pointer_real), dimension(:), allocatable, save :: send_nodes_values, recv_nodes_values
  integer, dimension(:), allocatable, save  :: request_send, request_recv
  integer, dimension(:,:), allocatable, save :: status_send, status_recv
  logical, dimension(:), allocatable, save :: partition_intersection_recv, partition_intersection_send
  integer, dimension(:), allocatable, save :: number_of_elements_to_receive, number_of_elements_to_send

  integer(kind = c_int8_t), dimension(:), pointer, save :: buffer_mpi

  ! Data
  integer(kind = c_int8_t), dimension(:), allocatable, save :: data
  integer, dimension(2), save                               :: data_b_header
  real, dimension(:), allocatable, save                     :: data_b_data
  integer(kind = c_int8_t), dimension(:), allocatable, save :: ldata

  integer, save :: mpi_comm, nprocs, rank

#include "mpif.h"

contains

                                                       ! TODO: Remove this unused argument
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
      subroutine donor_ele_data(eles, data)
        use iso_c_binding, only : c_int8_t
        implicit none
        integer, dimension(:), intent(in)                   :: eles
        integer(kind = c_int8_t), dimension(:), allocatable :: data
      end subroutine donor_ele_data

      subroutine unpack_data_b(data_b)
        use iso_c_binding, only : c_int8_t
        implicit none
        integer(kind = c_int8_t), dimension(:), intent(in) :: data_b
      end subroutine unpack_data_b
      
      subroutine intersection_calculation(positions_c, ele_a, ele_b, local)
        use iso_c_binding, only : c_int8_t
        implicit none
        ! dim x loc_c x nelements_c
        real, dimension(:, :, :), intent(in) :: positions_c
        integer, intent(in)     :: ele_a
        integer, intent(in)     :: ele_b
        logical, intent(in)     :: local
      end subroutine intersection_calculation
    end interface

    integer                     :: ierr, sends, recvs

    if(present(comm)) then
      mpi_comm = comm
    else
      mpi_comm = MPI_COMM_WORLD
    end if

    CALL MPI_Comm_rank(mpi_comm, rank, ierr)
    if(ierr /= MPI_SUCCESS) then
      FLAbort("MPI_Comm_rank error")
    end if
    CALL MPI_Comm_size(mpi_comm, nprocs, ierr)
    if(ierr /= MPI_SUCCESS) then
      FLAbort("MPI_Comm_size error")
    end if

    ! 0. Allocate and initialise arrays
    call initialise_parallel_supermesh(nprocs)

    ! 1. Calculate and communicate bounding box data for donor and target
    call step_1(positions_a, enlist_a, ele_owner_a, positions_b, enlist_b, ele_owner_b)

    ! 2. Use bounding box data to cull donor mesh
    call step_2(positions_b, enlist_b, ele_owner_b, sends, recvs, donor_ele_data)

    ! 5. Supermesh and call user specified element unpack and calculation functions
    call step_3(positions_a, enlist_a, ele_owner_a, positions_b, enlist_b, ele_owner_b, &
             &  sends, recvs, unpack_data_b, intersection_calculation)

    ! 6. Deallocate
    call finalise_parallel_supermesh()

  end subroutine parallel_supermesh


  subroutine step_1(positions_a, enlist_a, ele_owner_a, positions_b, enlist_b, ele_owner_b)
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

    integer  :: ierr

    allocate(bbox_a(2, size(positions_a,1)))
    bbox_a = partition_bbox(positions_a, enlist_a, ele_owner_a, rank)

    allocate(bbox_b(2, size(positions_b,1)))
    bbox_b = partition_bbox(positions_b, enlist_b, ele_owner_b, rank)

    allocate(parallel_bbox_a(2, size(positions_a,1), 0:nprocs-1))
    allocate(parallel_bbox_b(2, size(positions_b,1), 0:nprocs-1))

    ! Bounding box all-to-all
    call MPI_Allgather(bbox_a, 2 * size(positions_a,1), MPI_DOUBLE_PRECISION, &
      & parallel_bbox_a, 2 * size(positions_a,1), MPI_DOUBLE_PRECISION,    &
      & mpi_comm, ierr)
    if(ierr /= MPI_SUCCESS) then
      FLAbort("Unable to MPI_Allgather bbox_a(s) to parallel_bbox_a.")
    end if

    call MPI_Allgather(bbox_b, 2 * size(positions_b,1), MPI_DOUBLE_PRECISION, &
      & parallel_bbox_b, 2 * size(positions_b,1), MPI_DOUBLE_PRECISION,    &
      & mpi_comm, ierr)
    if(ierr /= MPI_SUCCESS) then
      FLAbort("Unable to MPI_Allgather bbox_b(s) to parallel_bbox_b.")
    end if

  end subroutine step_1


  subroutine step_2(positions_b, enlist_b, ele_owner_b, sends, recvs, donor_ele_data)
    ! dim x nnodes_b
    real, dimension(:, :), intent(in)    :: positions_b
    ! loc_b x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    ! elements_b
    integer, dimension(:), intent(in)    :: ele_owner_b
    integer, intent(out)                 :: sends
    integer, intent(out)                 :: recvs

    integer    :: ele_B, i, j, k, ierr, buffer_size, dp_extent, int_extent, position
    
    integer, dimension(:), allocatable :: elements_send
    integer :: nelements_send
    type(integer_set) :: nodes_send
    type(integer_hash_table) :: node_map
    integer :: nnodes_send
    integer :: node

    interface
      subroutine donor_ele_data(eles, data)
        use iso_c_binding, only : c_int8_t
        implicit none
        integer, dimension(:), intent(in)                   :: eles
        integer(kind = c_int8_t), dimension(:), allocatable :: data
      end subroutine donor_ele_data
    end interface

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

    allocate(send_nodes_connectivity(0:nprocs-1))
    do i=0,size(send_nodes_connectivity(:))-1
      nullify(send_nodes_connectivity(i)%p)
    end do

    allocate(recv_nodes_connectivity(0:nprocs-1))
    do i=0,size(recv_nodes_connectivity(:))-1
      nullify(recv_nodes_connectivity(i)%p)
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
    if(ierr /= MPI_SUCCESS) then
      FLAbort("MPI_TYPE_EXTENT double precision error")
    end if
    call MPI_TYPE_EXTENT(MPI_INTEGER, int_extent, ierr)
    if(ierr /= MPI_SUCCESS) then
      FLAbort("MPI_TYPE_EXTENT integer error")
    end if

    allocate(partition_intersection_recv(0:nprocs - 1), partition_intersection_send(0:nprocs - 1))
    partition_intersection_recv = .false.
    partition_intersection_send = .false.
    sends = -1
    recvs = -1
    do i = 0, nprocs - 1
      if(i == rank) cycle

      if(bboxes_intersect(bbox_a, parallel_bbox_b(:, :, i))) then
        ! Receiving data from process i
        partition_intersection_recv(i) = .true.
      end if

      if(bboxes_intersect(bbox_b, parallel_bbox_a(:, :, i))) then
        ! Sending data to process i
        partition_intersection_send(i) = .true.

        allocate(elements_send(size(enlist_b, 2)))
        nelements_send = 0
        do ele_B = 1, size(enlist_b, 2)
          ! Don't send non-owned elements
          if(ele_owner_b(ele_B) /= rank) cycle
          ! Don't send elements with non-intersecting bounding boxes
          if(.not. bboxes_intersect(bbox(positions_b(:, enlist_b(:, ele_B))), parallel_bbox_a(:, :, i))) cycle

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

        number_of_elements_and_nodes_to_send(1, i) = nelements_send
        number_of_elements_and_nodes_to_send(2, i) = nnodes_send

        allocate(send_element_uns(i)%p(nelements_send))
        send_element_uns(i)%p = elements_send(1:nelements_send)

        allocate(send_nodes_connectivity(i)%p(number_of_elements_and_nodes_to_send(1, i) * size(enlist_b, 1)))

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
            send_nodes_connectivity(i)%p((j - 1) * size(enlist_b, 1) + k) = fetch(node_map, enlist_b(k, elements_send(j)))
          end do
        end do
        call deallocate(node_map)
        deallocate(elements_send)

        allocate(send_nodes_values(i)%p(nnodes_send * size(positions_b, 1)))
#ifndef NDEBUG
        send_nodes_values(i)%p = -22.0D0
#endif
        do j = 1, nnodes_send
          node = fetch(nodes_send, j)
          send_nodes_values(i)%p((j - 1) * size(positions_b, 1) + 1:j * size(positions_b, 1)) = positions_b(:, node)
        end do
        call deallocate(nodes_send)

      ! ### Mesh send buffer allocation ###
        if ( number_of_elements_and_nodes_to_send(1,i) > 0 ) then
          call donor_ele_data(send_element_uns(i)%p, data)

          buffer_size = int_extent + int_extent +                          &
                    &   (size(send_nodes_connectivity(i)%p) * int_extent) + &
                    &   (size(send_nodes_values(i)%p)      * dp_extent ) + &
                    &   int_extent + (size(data))

          allocate(buffer_mpi(buffer_size))
          position = 0

          call MPI_Pack(number_of_elements_and_nodes_to_send(1, i), &
               & 1,                                                 &
               & MPI_INTEGER, buffer_mpi, buffer_size, position, MPI_COMM_WORLD, ierr)
          if(ierr /= MPI_SUCCESS) then
            FLAbort("MPI_Pack number of elements error")
          end if
          call MPI_Pack(number_of_elements_and_nodes_to_send(2, i), &
               & 1,                                                 &
               & MPI_INTEGER, buffer_mpi, buffer_size, position, MPI_COMM_WORLD, ierr)
          if(ierr /= MPI_SUCCESS) then
            FLAbort("MPI_Pack number of nodes error")
          end if
          call MPI_Pack(send_nodes_connectivity(i)%p,                &
               & size(send_nodes_connectivity(i)%p),                 &
               & MPI_INTEGER, buffer_mpi, buffer_size, position, MPI_COMM_WORLD, ierr)
          if(ierr /= MPI_SUCCESS) then
            FLAbort("MPI_Pack connectivity error")
          end if
          call MPI_Pack(send_nodes_values(i)%p,                     &
               & size(send_nodes_values(i)%p),                      &
               & MPI_DOUBLE_PRECISION, buffer_mpi, buffer_size, position, MPI_COMM_WORLD, ierr)
          if(ierr /= MPI_SUCCESS) then
            FLAbort("MPI_Pack values error")
          end if
          call MPI_Pack(size(data),                                 &
               & 1,                                                 &
               & MPI_INTEGER, buffer_mpi, buffer_size, position, MPI_COMM_WORLD, ierr)
          if(ierr /= MPI_SUCCESS) then
            FLAbort("MPI_Pack size of data error")
          end if
          call MPI_Pack(data,                                       &
               & size(data),                                        &
               & MPI_BYTE, buffer_mpi, buffer_size, position, MPI_COMM_WORLD, ierr)
          if(ierr /= MPI_SUCCESS) then
            FLAbort("MPI_Pack data error")
          end if

          allocate(send_buffer(i)%p(buffer_size))
          send_buffer(i)%p = buffer_mpi

          deallocate(data)
        else
          buffer_size = int_extent + int_extent

          allocate(buffer_mpi(buffer_size))
          position = 0
          call MPI_Pack(number_of_elements_and_nodes_to_send(1, i), &
               & 1,                                                 &
               & MPI_INTEGER, buffer_mpi, buffer_size, position, MPI_COMM_WORLD, ierr)
          if(ierr /= MPI_SUCCESS) then
            FLAbort("MPI_Pack number of elements error")
          end if
          call MPI_Pack(number_of_elements_and_nodes_to_send(2, i), &
               & 1,                                                 &
               & MPI_INTEGER, buffer_mpi, buffer_size, position, MPI_COMM_WORLD, ierr)
          if(ierr /= MPI_SUCCESS) then
            FLAbort("MPI_Pack number of nodes error")
          end if

          allocate(send_buffer(i)%p(buffer_size))
          send_buffer(i)%p = buffer_mpi
        end if

        sends = sends + 1

        call MPI_Isend(send_buffer(i)%p, &
               & size(send_buffer(i)%p), &
               & MPI_PACKED,             &
               & i, 0, MPI_COMM_WORLD, request_send(sends), ierr)
        if(ierr /= MPI_SUCCESS) then
          FLAbort("MPI_Isend error")
        end if
        deallocate(buffer_mpi)
    end if
  end do

  end subroutine step_2

  subroutine step_3(positions_a, enlist_a, ele_owner_a, &
                &   positions_b, enlist_b,       ele_owner_b, &
                &   sends, recvs, &
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
    integer, intent(inout)                  :: sends
    integer, intent(inout)                  :: recvs

    integer                                 :: i, j, k, l, m, n_C, ele_A, ele_B, nintersections, ierr, buffer_size, icount, position
    integer                                 :: status(MPI_STATUS_SIZE), nelements, nnodes
    real, dimension(:, :), allocatable      :: nodes_A, nodes_B
    real, dimension(:, :, :), allocatable   :: positions_c

    interface
      subroutine unpack_data_b(data_b)
        use iso_c_binding, only : c_int8_t
        implicit none
        integer(kind = c_int8_t), dimension(:), intent(in) :: data_b
      end subroutine unpack_data_b

      subroutine intersection_calculation(positions_c, ele_a, ele_b, local)
        use iso_c_binding, only : c_int8_t
        implicit none
        ! dim x loc_c x nelements_c
        real, dimension(:, :, :), intent(in) :: positions_c
        integer, intent(in) :: ele_a
        integer, intent(in) :: ele_b
        logical, intent(in) :: local
      end subroutine intersection_calculation
    end interface

    select case(size(positions_a, 1))
      case(1)
       ! 1D intervals
       allocate(positions_c(1, 2, 5), &
              & nodes_A(1, 2), nodes_B(1, 2))
      case(2)
        ! 2D triangles
        allocate(positions_c(2, 3, tri_buf_size), &
               & nodes_A(2, 3), nodes_B(2, 3))
      case(3)
        ! 3D tets
        allocate(positions_c(3, 4, tet_buf_size), &
               & nodes_A(3, 4), nodes_B(3, 4))
      case default
        FLAbort("Invalid dimension")
    end select
    nodes_A = 0.0
    nodes_B = 0.0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Parallel self-self runtime test !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

        call intersect_elements(nodes_A, nodes_B, n_C, positions_c)
        call intersection_calculation(positions_c(:, :, :n_C), ele_A, ele_B, .true.)

      end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Parallel self-other runtime test !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 0, nprocs - 1
      if(i == rank) cycle

      if(partition_intersection_recv(i)) then
        call MPI_Probe(i, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        if(ierr /= MPI_SUCCESS) then
          FLAbort("MPI_Probe error")
        end if
        call MPI_Get_Count(status, MPI_PACKED, icount, ierr)
        if(ierr /= MPI_SUCCESS) then
          FLAbort("MPI_Get_Count error")
        end if

        allocate(recv_buffer(i)%p(icount))
        call MPI_Recv (recv_buffer(i)%p, size(recv_buffer(i)%p), MPI_PACKED, & 
                  &  i, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        if(ierr /= MPI_SUCCESS) then
          FLAbort("MPI_Recv error")
        end if

        buffer_size = size(recv_buffer(i)%p)
        position = 0
        call MPI_UnPack ( recv_buffer(i)%p, buffer_size,  &
             & position, nelements,                        &
             & 1,                              &
             & MPI_INTEGER, MPI_COMM_WORLD, IERR)
        if(ierr /= MPI_SUCCESS) then
          FLAbort("MPI_UnPack number of elements error")
        end if
        allocate(recv_nodes_connectivity(i)%p(nelements * (size(positions_b,1) + 1)))

        call MPI_UnPack ( recv_buffer(i)%p, buffer_size,  &
             & position, nnodes,                    &
             & 1,                              &
             & MPI_INTEGER, MPI_COMM_WORLD, IERR)
        if(ierr /= MPI_SUCCESS) then
          FLAbort("MPI_UnPack number of nodes error")
        end if

        allocate(recv_nodes_values(i)%p(nnodes * (size(positions_b,1))))
        number_of_elements_and_nodes_to_receive(1, i) = nelements
        number_of_elements_and_nodes_to_receive(2, i) = nnodes

        if(number_of_elements_and_nodes_to_receive(1, i) <= 0) cycle

        call MPI_UnPack ( recv_buffer(i)%p, buffer_size,  &
             & position, recv_nodes_connectivity(i)%p,     &
             & size(recv_nodes_connectivity(i)%p),         &
             & MPI_INTEGER, MPI_COMM_WORLD, IERR)
        if(ierr /= MPI_SUCCESS) then
          FLAbort("MPI_UnPack connectivity error")
        end if

        call MPI_UnPack ( recv_buffer(i)%p, buffer_size,  &
             & position, recv_nodes_values(i)%p,          &
             & size(recv_nodes_values(i)%p),              &
             & MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
        if(ierr /= MPI_SUCCESS) then
          FLAbort("MPI_UnPack values error")
        end if

        call MPI_UnPack ( recv_buffer(i)%p, buffer_size,  &
             & position, k,                               &
             & 1,                                         &
             & MPI_INTEGER, MPI_COMM_WORLD, IERR)
        if(ierr /= MPI_SUCCESS) then
          FLAbort("MPI_UnPack size of data error")
        end if

        allocate(ldata(k))
        call MPI_UnPack ( recv_buffer(i)%p, buffer_size,  &
             & position, ldata,                           &
             & k,                                         &
             & MPI_BYTE, MPI_COMM_WORLD, IERR)
        if(ierr /= MPI_SUCCESS) then
          FLAbort("MPI_UnPack data error")
        end if

        call unpack_data_b(ldata)

        do ele_B = 1, nelements
          do k = 1, size(enlist_b, 1)
            do j = 1, size(positions_b, 1)
              l = recv_nodes_connectivity(i)%p((ele_B - 1) * size(enlist_b, 1) + k)            
              nodes_B(j, k) = recv_nodes_values(i)%p((l - 1) * size(positions_b, 1) + j)
            end do
          end do

          call rtree_intersection_finder_find(nodes_B)
          call rtree_intersection_finder_query_output(nintersections)
          do k = 1, nintersections
            call rtree_intersection_finder_get_output(ele_A, k)
            if(ele_owner_a(ele_A) /= rank) cycle
            nodes_A = positions_a(:, enlist_a(:, ele_A))

            call intersect_elements(nodes_A, nodes_B, n_C, positions_c)
            call intersection_calculation(positions_c(:, :, :n_C), ele_A, ele_B, .false.)
          end do
        end do

        deallocate(recv_nodes_values(i)%p)
        deallocate(ldata)
      end if
    end do

    deallocate(positions_c, nodes_A, nodes_B)

    sends = sends + 1
    call MPI_Waitall(sends, request_send(0:sends), status_send(:,0:sends), ierr)
    if(ierr /= MPI_SUCCESS) then
      FLAbort("Unable to setup MPI_Waitall(send).")
    end if

  end subroutine step_3

  subroutine initialise_parallel_supermesh(nprocs)
    integer, intent(in)       :: nprocs

    integer  :: i

    parallel_supermesh_allocated = .true.
    allocate(number_of_elements_to_receive(0:nprocs - 1), &
           & number_of_elements_to_send(0:nprocs -1),     &
           & request_send(0: nprocs - 1), request_recv(0: nprocs - 1), &
           & status_send(MPI_STATUS_SIZE, 0: nprocs - 1), status_recv(MPI_STATUS_SIZE, 0: nprocs - 1))
    number_of_elements_to_receive = -11
    number_of_elements_to_send    = -11
    forall(i = 0:nprocs - 1)
      status_send(:, i) = MPI_STATUS_IGNORE
      status_recv(:, i) = MPI_STATUS_IGNORE
    end forall
    request_send = MPI_REQUEST_NULL
    request_recv = MPI_REQUEST_NULL

    allocate(number_of_elements_and_nodes_to_receive(2,0:nprocs - 1), number_of_elements_and_nodes_to_send(2,0:nprocs - 1))
    number_of_elements_and_nodes_to_receive = 0
    number_of_elements_and_nodes_to_send = 0
  end subroutine initialise_parallel_supermesh

  subroutine finalise_parallel_supermesh()
    integer :: i, ntests

    if(parallel_supermesh_allocated) then
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

      do i=0,size(send_nodes_connectivity(:))-1
        if ( .NOT. ASSOCIATED(send_nodes_connectivity(i)%p) ) cycle
        deallocate(send_nodes_connectivity(i)%p)
      end do
      deallocate(send_nodes_connectivity)

      do i=0,size(recv_nodes_connectivity(:))-1
        if ( .NOT. ASSOCIATED(recv_nodes_connectivity(i)%p) ) cycle
        deallocate(recv_nodes_connectivity(i)%p)
      end do
      deallocate(recv_nodes_connectivity)

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

      deallocate(parallel_bbox_a, parallel_bbox_b, bbox_a, bbox_b)
      deallocate(request_send, request_recv, status_send, status_recv, partition_intersection_send, partition_intersection_recv)
      deallocate(number_of_elements_to_receive, number_of_elements_to_send)
      deallocate(number_of_elements_and_nodes_to_receive, number_of_elements_and_nodes_to_send)

      parallel_supermesh_allocated = .false.
    end if
    call cintersection_finder_reset(ntests)

  end subroutine finalise_parallel_supermesh

end module libsupermesh_parallel_supermesh
