#include "fdebug.h"

module libsupermesh_halo_ownership

  use libsupermesh_debug
  use libsupermesh_read_halos

  implicit none

#include "mpif.h"

  private

  public :: element_ownership, node_ownership, universal_node_numbering

contains

  subroutine element_ownership(nnodes, enlist, halo, ele_owner)
    integer, intent(in) :: nnodes
    ! loc x nelements
    integer, dimension(:, :), intent(in) :: enlist
    ! level 2 halo
    type(halo_type), intent(in) :: halo
    ! nelements
    integer, dimension(:), intent(out) :: ele_owner

    integer :: i
    integer, dimension(:), allocatable :: node_owner

    allocate(node_owner(nnodes))
    call node_ownership(nnodes, halo, node_owner)
    
    do i = 1, size(enlist, 2)
      ele_owner(i) = minval(node_owner(enlist(:, i)))
    end do

    deallocate(node_owner)

  end subroutine element_ownership

  subroutine node_ownership(nnodes, halo, node_owner)
    integer, intent(in) :: nnodes
    ! level 2 halo
    type(halo_type), intent(in) :: halo
    ! nelements
    integer, dimension(:), intent(out) :: node_owner

    integer :: i, ierr, process

    call MPI_Comm_rank(MPI_COMM_WORLD, process, ierr)
    assert(ierr == MPI_SUCCESS)
    node_owner(:halo%npnodes) = process
    node_owner(halo%npnodes + 1:) = huge(0)

    do i = 1, halo%nprocs
      if(size(halo%recv(i)%val) > 0) then
        assert(all(halo%recv(i)%val > halo%npnodes))
        assert(all(halo%recv(i)%val <= nnodes))
        node_owner(halo%recv(i)%val) = i - 1
      end if
    end do

  end subroutine node_ownership

  subroutine universal_node_numbering(halo, uns)
    ! level 2 halo
    type(halo_type), intent(in) :: halo
    ! nnodes
    integer, dimension(:), intent(out) :: uns

    integer :: i, ierr, un_offset
    integer, dimension(:), allocatable :: recv_types, send_types, requests, statuses

    call MPI_Scan(halo%npnodes, un_offset, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    assert(ierr == MPI_SUCCESS)
    un_offset = un_offset - halo%npnodes

    do i = 1, halo%npnodes
      uns(i) = i + un_offset
    end do

    allocate(send_types(halo%nprocs), recv_types(halo%nprocs), requests(2 * halo%nprocs))
    send_types = MPI_DATATYPE_NULL
    recv_types = MPI_DATATYPE_NULL
    requests = MPI_REQUEST_NULL
    do i = 1, halo%nprocs      
      if(size(halo%send(i)%val) > 0) then
        call MPI_Type_create_indexed_block(size(halo%send(i)%val), 1, halo%send(i)%val - 1, MPI_INTEGER, send_types(i), ierr)
        assert(ierr == MPI_SUCCESS)
        call MPI_Type_commit(send_types(i), ierr)
        assert(ierr == MPI_SUCCESS)

        call MPI_Isend(uns, 1, send_types(i), i - 1, halo%process, MPI_COMM_WORLD, requests(i), ierr)
        assert(ierr == MPI_SUCCESS)
      end if

      if(size(halo%recv(i)%val) > 0) then
        call MPI_Type_create_indexed_block(size(halo%recv(i)%val), 1, halo%recv(i)%val - 1, MPI_INTEGER, recv_types(i), ierr)
        assert(ierr == MPI_SUCCESS)
        call MPI_Type_commit(recv_types(i), ierr)
        assert(ierr == MPI_SUCCESS)

        call MPI_Irecv(uns, 1, recv_types(i), i - 1, i - 1, MPI_COMM_WORLD, requests(halo%nprocs + i), ierr)
        assert(ierr == MPI_SUCCESS)
      end if
    end do

    allocate(statuses(MPI_STATUS_SIZE * size(requests)))
    call MPI_Waitall(size(requests), requests, statuses, ierr)
    deallocate(requests, statuses)
    assert(ierr == MPI_SUCCESS)

    do i = 1, halo%nprocs
      if(send_types(i) /= MPI_DATATYPE_NULL) then
        call MPI_Type_free(send_types(i), ierr)
        assert(ierr == MPI_SUCCESS)
      end if
      if(recv_types(i) /= MPI_DATATYPE_NULL) then
        call MPI_Type_free(recv_types(i), ierr)
        assert(ierr == MPI_SUCCESS)
      end if
    end do
    deallocate(send_types, recv_types)

  end subroutine universal_node_numbering

end module libsupermesh_halo_ownership
