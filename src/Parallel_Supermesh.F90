#include "fdebug.h"

module libsupermesh_parallel_supermesh

  use iso_fortran_env, only : output_unit
  use libsupermesh_integer_hash_table
  use libsupermesh_fldebug
  use libsupermesh_read_halos
  use libsupermesh_intersection_finder
  use libsupermesh_tri_intersection_module
  use libsupermesh_tet_intersection_module
  use libsupermesh_construction
  use libsupermesh_fields, only : triangle_area

  implicit none

  private

  public :: parallel_supermesh

  type pointer_real
    real, dimension(:), pointer :: p
  end type pointer_real

  type pointer_integer
    integer, dimension(:), pointer :: p
  end type pointer_integer

  logical, save :: parallel_supermesh_allocated = .false.
  real, dimension(:,:), allocatable, save :: bbox_a, bbox_b
  real, dimension(:,:,:), allocatable, save :: parallel_bbox_a, parallel_bbox_b

  type(pointer_integer), dimension(:), allocatable, save :: send_element_uns

  type(pointer_real), dimension(:), allocatable, save :: recv_buffer, send_buffer
  integer, dimension(:), allocatable, save  :: request_send, request_recv
  integer, dimension(:,:), allocatable, save :: status_send, status_recv
  logical, dimension(:), allocatable, save :: partition_intersection_recv, partition_intersection_send
  integer, dimension(:), allocatable, save :: number_of_elements_to_receive, number_of_elements_to_send, temp_elements_uns

  integer :: mpi_comm, nprocs, rank

#include "mpif.h"

contains

  subroutine parallel_supermesh(positions_a, enlist_a, un_a, ele_owner_a, &
                             &  positions_b, enlist_b,       ele_owner_b, &
                             &  donor_ele_data, intersection_calculation, &
                             &  comm)
    ! dim x nnodes_a
    real, dimension(:, :), intent(in)    :: positions_a
    ! loc_a x nelements_a
    integer, dimension(:, :), intent(in) :: enlist_a
    ! nnodes_a
    integer, dimension(:), intent(in)    :: un_a
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
      subroutine donor_ele_data(eles, data, ndata)
        use iso_c_binding, only : c_ptr
        implicit none
        integer, dimension(:), intent(in) :: eles
        type(c_ptr), intent(out)          :: data
        integer, intent(out)              :: ndata
      end subroutine donor_ele_data

      subroutine intersection_calculation(positions_c, ele_a, un_a, n_C, ele_b, data_a, ndata_a)
        use iso_c_binding, only : c_ptr
        implicit none
        ! dim x loc_c x nelements_c
        real, dimension(:, :, :), intent(in) :: positions_c
        integer, intent(in)     :: ele_a
        integer, intent(in)     :: un_a
        integer, intent(in)     :: n_C
        integer, intent(in)     :: ele_b
        type(c_ptr), intent(in) :: data_a
        integer, intent(in)     :: ndata_a
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
    call step_2(positions_b, enlist_b, ele_owner_b, sends, recvs)

    ! 3. Pack non-culled mesh data for communication, calling user specified data functions
    call step_3(positions_b, enlist_b, sends)

    ! 4. Communicate donor mesh and mesh data
    call step_4(positions_b, sends, recvs)

    ! 5. Supermesh and call user specified element unpack and calculation functions
    call step_5(positions_a, enlist_a, un_a, ele_owner_a, positions_b, enlist_b, ele_owner_b, sends, recvs, intersection_calculation)

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

    allocate(bbox_a(size(positions_a,1), 2))
    bbox_a = partition_bbox(positions_a, enlist_a, ele_owner_a, rank)

    allocate(bbox_b(size(positions_b,1), 2))
    bbox_b = partition_bbox(positions_b, enlist_b, ele_owner_b, rank)

    allocate(parallel_bbox_a(2, size(positions_a,1), 0:nprocs-1))
    allocate(parallel_bbox_b(2, size(positions_b,1), 0:nprocs-1))

    ! Bounding box all-to-all
    call MPI_Allgather(bbox_a, 4, MPI_DOUBLE_PRECISION, &
      & parallel_bbox_a, 4, MPI_DOUBLE_PRECISION,    &
      & mpi_comm, ierr)
    if(ierr /= MPI_SUCCESS) then
      FLAbort("Unable to MPI_Allgather bbox_a(s) to parallel_bbox_a.")
    end if

    call MPI_Allgather(bbox_b, 4, MPI_DOUBLE_PRECISION, &
      & parallel_bbox_b, 4, MPI_DOUBLE_PRECISION,    &
      & mpi_comm, ierr)
    if(ierr /= MPI_SUCCESS) then
      FLAbort("Unable to MPI_Allgather bbox_b(s) to parallel_bbox_b.")
    end if

  end subroutine step_1


  subroutine step_2(positions_b, enlist_b, ele_owner_b, sends, recvs)
    ! dim x nnodes_b
    real, dimension(:, :), intent(in)    :: positions_b
    ! loc_b x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    ! elements_b
    integer, dimension(:), intent(in)    :: ele_owner_b
    integer, intent(out)                 :: sends
    integer, intent(out)                 :: recvs

    integer    :: ele_B, i, l, ierr

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

          call MPI_Irecv(number_of_elements_to_receive(i), 1, MPI_INTEGER, &
              & i, i, mpi_comm, &
              & request_recv(recvs), ierr)
          if(ierr /= MPI_SUCCESS) then
            FLAbort("Unable to setup MPI_Irecv.")
          end if
        end if
      end if

      if(bboxes_intersect(bbox_b, parallel_bbox_a(:, :, i))) then
        if(i /= rank) then
          sends = sends + 1
          partition_intersection_send(i) = .true.
          allocate(temp_elements_uns(size(enlist_b, 2)))
          do ele_B = 1, size(enlist_b, 2)
            if(ele_owner_b(ele_B) /= rank) cycle
            if(.not. bboxes_intersect(bbox(positions_b(:, enlist_b(:, ele_B))), parallel_bbox_a(:, :, i))) cycle
            l = l + 1                     ! Keep a counter
            temp_elements_uns(l) = ele_B  ! Keep the actual element
          end do
          number_of_elements_to_send(i) = l
          allocate(send_element_uns(i)%p(l))
          send_element_uns(i)%p = temp_elements_uns(1:l)
          deallocate(temp_elements_uns)
          call MPI_Isend(number_of_elements_to_send(i), 1, MPI_INTEGER, i, rank, &
              & mpi_comm, request_send(sends), ierr)
          if(ierr /= MPI_SUCCESS) then
            FLAbort("Unable to setup MPI_Isend.")
          end if

          ! ### Mesh send buffer allocation ###
          allocate(send_buffer(i)%p(l * (size(positions_b,1) + 1) * size(positions_b,1)))
        end if
      end if
    end do

    recvs = recvs + 1
    call MPI_Waitall(recvs, request_recv(0:recvs), status_recv(:,0:recvs), ierr)
    if(ierr /= MPI_SUCCESS) then
      FLAbort("Unable to setup MPI_Waitall(recv).")
    end if

    do i = 0, nprocs - 1
      if(partition_intersection_recv(i)) then
        ! ### Mesh receivebuffer allocation ###
        allocate(recv_buffer(i)%p( number_of_elements_to_receive(i) * (size(positions_b,1) + 1) * size(positions_b,1)))
      end if
    end do

  end subroutine step_2


  subroutine step_3(positions_b, enlist_b, sends)
    ! dim x nnodes_b
    real, dimension(:, :), intent(in)    :: positions_b
    ! loc_b x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    integer, intent(inout)               :: sends

    integer                              :: j, k, l, ierr
    type(tri_type)                       :: tri_B

    do j = 0, nprocs - 1
      if(.not. associated(send_element_uns(j)%p)) cycle
      if(size(send_element_uns(j)%p) == 0) cycle
      k = 1

      do l = 1, size(send_element_uns(j)%p)
        tri_B%v = positions_b(:, enlist_b(:, send_element_uns(j)%p(l)))
        send_buffer(j)%p(k) = tri_B%v(1, 1)
        k = k + 1
        send_buffer(j)%p(k) = tri_B%v(2, 1)
        k = k + 1
        send_buffer(j)%p(k) = tri_B%v(1, 2)
        k = k + 1
        send_buffer(j)%p(k) = tri_B%v(2, 2)
        k = k + 1
        send_buffer(j)%p(k) = tri_B%v(1, 3)
        k = k + 1
        send_buffer(j)%p(k) = tri_B%v(2, 3)
        k = k + 1
      end do
    end do

    sends = sends + 1
    call MPI_Waitall(sends, request_send(0:sends), status_send(:,0:sends), ierr)
    if(ierr /= MPI_SUCCESS) then
      FLAbort("Unable to setup MPI_Waitall(send).")
    end if

  end subroutine step_3


  subroutine step_4(positions_b, sends, recvs)
    ! dim x nnodes_b
    real, dimension(:, :), intent(in)  :: positions_b
    integer, intent(inout)             :: sends
    integer, intent(inout)             :: recvs

    integer                            :: i, j, ierr

    forall(i = 0: nprocs - 1)
      status_send(:, i) = MPI_STATUS_IGNORE
      status_recv(:, i) = MPI_STATUS_IGNORE
    end forall
    request_send = MPI_REQUEST_NULL
    request_recv = MPI_REQUEST_NULL

    sends = -1
    recvs = -1
    do j = 0,nprocs - 1
      if(rank == j)cycle

      if(partition_intersection_recv(j) .and. (number_of_elements_to_receive(j) /=0) ) then
        recvs = recvs + 1

        call MPI_Irecv(recv_buffer(j)%p, &
          & number_of_elements_to_receive(j) * (size(positions_b,1) + 1) * size(positions_b,1), &
          & MPI_DOUBLE_PRECISION, &
          & j, MPI_ANY_TAG, mpi_comm, &
          & request_recv(recvs), ierr)
        if(ierr /= MPI_SUCCESS) then
          FLAbort("Unable to setup MPI_Irecv(data).")
        end if
      end if

      if( (partition_intersection_send(j)) .and. (size(send_buffer(j)%p) /=0) ) then
        sends = sends + 1

        call MPI_Isend(send_buffer(j)%p, &
          & size(send_buffer(j)%p), &
          & MPI_DOUBLE_PRECISION, &
          & j, 0, mpi_comm, request_send(sends), ierr)
        if(ierr /= MPI_SUCCESS) then
          FLAbort("Unable to setup MPI_Isend(data).")
        end if
      end if
    end do

  end subroutine step_4


  subroutine step_5(positions_a, enlist_a, un_a, ele_owner_a, &
                &   positions_b, enlist_b,       ele_owner_b, &
                &   sends, recvs, &
                &   intersection_calculation)
    ! dim x nnodes_a
    real, dimension(:, :), intent(in)       :: positions_a
    ! loc_a x nelements_a
    integer, dimension(:, :), intent(in)    :: enlist_a
    ! nnodes_a
    integer, dimension(:), intent(in)       :: un_a
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

    integer                                 :: i, j, k, l, m, n_C, ele_A, ele_B, nintersections, ntests, ierr
    real, dimension(:, :), allocatable      :: nodes_A, nodes_B
    integer, dimension(:, :), allocatable   :: comm_enlist_B
    real, dimension(:, :), allocatable      :: comm_coords_B
    real, dimension(:, :, :), allocatable   :: positions_c

    ! JRM: FIXME
    type(c_ptr) :: dummy_data_a

    interface
      subroutine intersection_calculation(positions_c, ele_a, un_a, n_C, ele_b, data_a, ndata_a)
        use iso_c_binding, only : c_ptr
        implicit none
        ! dim x loc_c x nelements_c
        real, dimension(:, :, :), intent(in) :: positions_c
        integer, intent(in)     :: ele_a
        integer, intent(in)     :: un_a
        integer, intent(in)     :: n_C
        integer, intent(in)     :: ele_b
        type(c_ptr), intent(in) :: data_a
        integer, intent(in)     :: ndata_a
      end subroutine intersection_calculation
    end interface

     if (size(positions_a,1) .eq. 1) then
       ! 1D Vectors
       allocate(positions_c(1, 2, 2))
       allocate(nodes_A(1, 5), nodes_B(1, 2))
     else if (size(positions_a,1) .eq. 2) then
       if (size(enlist_a,1) .eq. 3) then
         ! 2D Triangles
         allocate(positions_c(2, 3, tri_buf_size))
         allocate(nodes_A(2, 3), nodes_B(2, 3))
       else if (size(enlist_a,1) .eq. 4) then
         ! 2D Tets
         allocate(positions_c(3, 4, tet_buf_size))
         allocate(nodes_A(3, 4), nodes_B(3, 4))
       end if
     else if (size(positions_a,1) .eq. 1) then
       ! 3D
     end if
     nodes_A = 0.0
     nodes_B = 0.0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Parallel self-self runtime test !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call rtree_intersection_finder_set_input(positions_b, enlist_b)
    do ele_A = 1, size(enlist_a, 2)
      if(ele_owner_a(ele_A) /= rank) cycle
      nodes_A = positions_a(:, enlist_a(:, ele_A))
      call rtree_intersection_finder_find(nodes_A)
      call rtree_intersection_finder_query_output(nintersections)
      do i = 1, nintersections
        call rtree_intersection_finder_get_output(ele_B, i)
        if(ele_owner_b(ele_B) /= rank) cycle
        nodes_B = positions_b(:, enlist_b(:, ele_B))

        call intersect_elements(nodes_A, nodes_B, n_C, positions_c)
        call intersection_calculation(positions_c, ele_A, un_a(ele_A), n_C, ele_B, dummy_data_a, 0)

      end do
    end do
    call rtree_intersection_finder_reset(ntests)

    sends = sends + 1
    recvs = recvs + 1
    call MPI_Waitall(sends, request_send(0:sends), status_send(:,0:sends), ierr)
    if(ierr /= MPI_SUCCESS) then
      FLAbort("Unable to setup MPI_Waitall(send).")
    end if
    call MPI_Waitall(recvs, request_recv(0:recvs), status_recv(:,0:recvs), ierr)
    if(ierr /= MPI_SUCCESS) then
      FLAbort("Unable to setup MPI_Waitall(recv).")
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Parallel self-other runtime test !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 0, nprocs - 1
      if(i == rank) cycle
      if(number_of_elements_to_receive(i) <= 0) cycle
      allocate(comm_coords_B( size(positions_b,1),      number_of_elements_to_receive(i) * (size(positions_b,1) + 1)), &
             & comm_enlist_B((size(positions_b,1) + 1), number_of_elements_to_receive(i)))
      m = 1
      do j = 1, number_of_elements_to_receive(i)
        do k = 1, (size(positions_b,1) + 1)
          do l = 1, size(positions_b,1)
            comm_coords_B(l, (j - 1) * (size(positions_b,1) + 1) + k) = recv_buffer(i)%p(size(positions_b,1) * (m - 1) + l)
          end do
          comm_enlist_B(k, j) = m
          m = m + 1
        end do
      end do
      call rtree_intersection_finder_set_input(comm_coords_B, comm_enlist_B)
      do ele_A = 1, size(enlist_a, 2)
        if(ele_owner_a(ele_A) /= rank) cycle
        nodes_A = positions_a(:, enlist_a(:, ele_A))
        call rtree_intersection_finder_find(nodes_A)
        call rtree_intersection_finder_query_output(nintersections)
        do k = 1, nintersections
          call rtree_intersection_finder_get_output(ele_B, k)
          nodes_B = comm_coords_B(:, (ele_B - 1) * (size(positions_b,1) + 1) + 1:ele_B * (size(positions_b,1) + 1))

          call intersect_elements(nodes_A, nodes_B, n_C, positions_c)
          call intersection_calculation(positions_c, ele_A, un_a(ele_A), n_C, ele_B, dummy_data_a, 0)

        end do
      end do
      call rtree_intersection_finder_reset(ntests)
      deallocate(comm_coords_B, &
               & comm_enlist_B)
    end do

    deallocate(positions_c, nodes_A, nodes_B)

  end subroutine step_5

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

      deallocate(parallel_bbox_a, parallel_bbox_b, bbox_a, bbox_b)
      deallocate(request_send, request_recv, status_send, status_recv, partition_intersection_send, partition_intersection_recv)
      deallocate(number_of_elements_to_receive, number_of_elements_to_send)

      parallel_supermesh_allocated = .false.
    end if
    call cintersection_finder_reset(ntests)

  end subroutine finalise_parallel_supermesh

end module libsupermesh_parallel_supermesh
