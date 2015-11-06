#include "fdebug.h"

module libsupermesh_parallel_supermesh

  use iso_fortran_env, only : output_unit
  use libsupermesh_integer_hash_table
  use libsupermesh_fldebug
  use libsupermesh_read_halos
  use libsupermesh_intersection_finder
  use libsupermesh_tri_intersection_module
  use libsupermesh_fields, only : triangle_area

  implicit none

  private

  public :: parallel_supermesh, finalise_parallel_supermesh

  type pointer_real
    real, dimension(:), pointer :: p
  end type pointer_real

  type pointer_integer
    integer, dimension(:), pointer :: p
  end type pointer_integer

  logical, save, private :: parallel_supermesh_allocated = .false.
  real, dimension(:,:), allocatable    :: bbox_a, bbox_b
  real, dimension(:,:,:), allocatable  :: parallel_bbox_a, parallel_bbox_b

  type(pointer_integer), dimension(:), allocatable  :: send_element_uns

  type(pointer_real), dimension(:), allocatable  :: recv_buffer, send_buffer
  integer, dimension(:), allocatable   :: request_send, request_recv
  integer, dimension(:,:), allocatable :: status_send, status_recv
  logical, dimension(:), allocatable   :: partition_intersection_recv, partition_intersection_send
  integer, dimension(:), allocatable   :: number_of_elements_to_receive, number_of_elements_to_send, temp_elements_uns

  integer                              :: nprocs, rank

#include "mpif.h"

contains

  subroutine parallel_supermesh(positions_a, enlist_a, un_a, ele_owner_a, &
                             &  positions_b, enlist_b,       ele_owner_b, &
                             &  donor_ele_data, unpack_donor_ele_data, intersection_calculation)
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
    interface
      subroutine donor_ele_data(eles, data, ndata)
        use iso_c_binding, only : c_ptr
        implicit none
        integer, dimension(:), intent(in) :: eles
        type(c_ptr), intent(out)          :: data
        integer, intent(out)              :: ndata
      end subroutine donor_ele_data

      subroutine unpack_donor_ele_data(ele, proc, data, ndata, ele_data, nele_data)
        use iso_c_binding, only : c_ptr
        implicit none
        integer, intent(in)      :: ele
        integer, intent(in)      :: proc
        type(c_ptr), intent(in)  :: data
        integer, intent(in)      :: ndata
        type(c_ptr), intent(out) :: ele_data
        integer, intent(out)     :: nele_data
      end subroutine unpack_donor_ele_data

      subroutine intersection_calculation(positions_c, ele_a, proc_a, ele_b, ele_data_a, nele_data_a)
        use iso_c_binding, only : c_ptr
        implicit none
        ! dim x loc_c x nelements_c
        real, dimension(:, :, :), intent(in) :: positions_c
        integer, intent(in)     :: ele_a
        integer, intent(in)     :: proc_a
        integer, intent(in)     :: ele_b
        type(c_ptr), intent(in) :: ele_data_a
        integer, intent(in)     :: nele_data_a
      end subroutine intersection_calculation
    end interface

    integer                     :: i, ierr, sends, recvs

! find out MY process ID, and how many processes were started.
    CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    if(ierr /= MPI_SUCCESS) then
      FLAbort("Unable to MPI_Comm_rank.")
    end if
    CALL MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
    if(ierr /= MPI_SUCCESS) then
      FLAbort("Unable to MPI_Comm_size.")
    end if

    ! 0. Allocate and initialise arrays
    parallel_supermesh_allocated = .true.
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

    ! 1. Calculate and Communicate bounding box data for donor and target
    !
    call step_1(positions_a, enlist_a, ele_owner_a, positions_b, enlist_b, ele_owner_b)

    ! 2. Use bounding box data to cull donor mesh
    !
    call step_2(positions_b, enlist_b, ele_owner_b, sends, recvs)

    ! 3. Pack non-culled mesh data for communication, calling user specified data functions
    !
    call step_3(positions_b, enlist_b, sends)

    ! 4. Communicate donor mesh and mesh data
    !
    call step_4(positions_b, sends, recvs)

    ! 5. Supermesh and call user specified element unpack and calculation functions
    !
    call step_5(positions_a, enlist_a, ele_owner_a, positions_b, enlist_b, ele_owner_b, sends, recvs)

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

    allocate(bbox_a(size(positions_a,1),2))
    bbox_a = partition_bbox_real(positions_a, enlist_a, ele_owner_a, rank)

    allocate(bbox_b(size(positions_b,1),2))
    bbox_b = partition_bbox_real(positions_b, enlist_b, ele_owner_b, rank)

    allocate(parallel_bbox_a(2, size(positions_a,1), 0:nprocs-1))
    allocate(parallel_bbox_b(2, size(positions_b,1), 0:nprocs-1))

    ! Bounding box all-to-all
    call MPI_Allgather(bbox_a, 4, MPI_DOUBLE_PRECISION, &
      & parallel_bbox_a, 4, MPI_DOUBLE_PRECISION,    &
      & MPI_COMM_WORLD, ierr)
    if(ierr /= MPI_SUCCESS) then
      FLAbort("Unable to MPI_Allgather bbox_a(s) to parallel_bbox_a.")
    end if

    call MPI_Allgather(bbox_b, 4, MPI_DOUBLE_PRECISION, &
      & parallel_bbox_b, 4, MPI_DOUBLE_PRECISION,    &
      & MPI_COMM_WORLD, ierr)
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
              & i, i, MPI_COMM_WORLD, &
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
              & MPI_COMM_WORLD, request_send(sends), ierr)
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
          & j, MPI_ANY_TAG, MPI_COMM_WORLD, &
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
          & j, 0, MPI_COMM_WORLD, request_send(sends), ierr)
        if(ierr /= MPI_SUCCESS) then
          FLAbort("Unable to setup MPI_Isend(data).")
        end if
      end if
    end do

  end subroutine step_4


  subroutine step_5(positions_a, enlist_a, ele_owner_a, positions_b, enlist_b, ele_owner_b, sends, recvs)
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

    integer                                 :: i, j, k, l, m, n_trisc, ele_A, ele_B, ele_C, nintersections, ntests, ierr
    real                                    :: area_parallel
    type(tri_type)                          :: tri_A, tri_B
    integer, dimension(:, :), allocatable   :: comm_enlist_B
    real, dimension(:, :), allocatable      :: comm_coords_B
    type(tri_type), dimension(tri_buf_size) :: trisC
#ifdef NDEBUG
    real, dimension(:), allocatable         :: areas_parallel
#endif

    area_parallel = 0.0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Parallel self-self runtime test !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call rtree_intersection_finder_set_input(positions_b, enlist_b)
    do ele_A = 1, size(enlist_a, 2)
      if(ele_owner_a(ele_A) /= rank) cycle
      tri_A%v = positions_a(:, enlist_a(:, ele_A))
      call rtree_intersection_finder_find(tri_A%v)
      call rtree_intersection_finder_query_output(nintersections)
      do i = 1, nintersections
        call rtree_intersection_finder_get_output(ele_B, i)
        if(ele_owner_b(ele_B) /= rank) cycle
        tri_B%v = positions_b(:, enlist_b(:, ele_B))

!        local_iter = local_iter + 1
        call intersect_tris(tri_A, tri_B, trisC, n_trisC)

        do ele_C = 1, n_trisC
          area_parallel = area_parallel + triangle_area(trisC(ele_C)%v)
!          local_iter_actual = local_iter_actual + 1
        end do
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
        tri_A%v = positions_a(:, enlist_a(:, ele_A))
        call rtree_intersection_finder_find(tri_A%v)
        call rtree_intersection_finder_query_output(nintersections)
        do k = 1, nintersections
          call rtree_intersection_finder_get_output(ele_B, k)
          tri_B%v = comm_coords_B(:, (ele_B - 1) * (size(positions_b,1) + 1) + 1:ele_B * (size(positions_b,1) + 1))

!          local_iter = local_iter + 1

          call intersect_tris(tri_A, tri_B, trisC, n_trisC)

          do ele_C = 1, n_trisC
            area_parallel = area_parallel + triangle_area(trisC(ele_C)%v)
!            local_iter_actual = local_iter_actual + 1
          end do
        end do
      end do
      call rtree_intersection_finder_reset(ntests)
      deallocate(comm_coords_B, &
               & comm_enlist_B)
    end do

#ifdef NDEBUG
    allocate(areas_parallel(0:nprocs - 1))
    write(output_unit, *) rank, ": area_parallel:",area_parallel
    call MPI_Gather(area_parallel, 1, MPI_DOUBLE_PRECISION, &
       & areas_parallel, 1, MPI_DOUBLE_PRECISION, &
       & 0, MPI_COMM_WORLD, ierr)
    if(ierr /= MPI_SUCCESS) then
      FLAbort("Unable to setup MPI_Gather(area_parallel).")
    end if

    if(rank == 0) then
      write(output_unit, "(i5,a,F19.15,a)") rank, ": Total parallel intersection area   : ", sum(areas_parallel)," ."
    end if

    deallocate(areas_parallel)
#endif

  end subroutine step_5

  subroutine finalise_parallel_supermesh()
    integer :: i

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

  end subroutine finalise_parallel_supermesh

end module libsupermesh_parallel_supermesh
