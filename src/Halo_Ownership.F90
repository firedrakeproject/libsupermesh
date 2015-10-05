#include "fdebug.h"

module libsupermesh_halo_ownership

  use libsupermesh_fldebug
  use libsupermesh_read_halos

  implicit none

#include "mpif.h"

  private

  public :: element_ownership, node_ownership

contains

  subroutine element_ownership(nnodes, enlist, halo, ele_owner)
    integer, intent(in) :: nnodes
    ! loc x nelements
    integer, dimension(:, :), intent(in) :: enlist
    ! level 2 halo
    type(halo_type), intent(in) :: halo
    ! nelements
    integer, dimension(:), intent(out) :: ele_owner

    integer :: i, ierr, process
    integer, dimension(:), allocatable :: node_owner

    call MPI_Comm_rank(MPI_COMM_WORLD, process, ierr)
    assert(ierr == MPI_SUCCESS)
    allocate(node_owner(nnodes))
    node_owner(:halo%npnodes) = process
    node_owner(halo%npnodes + 1:) = huge(0)

    do i = 1, halo%nprocs
      if(size(halo%recv(i)%val) > 0) then
        assert(all(halo%recv(i)%val > halo%npnodes))
        assert(all(halo%recv(i)%val <= nnodes))
        node_owner(halo%recv(i)%val) = i - 1
      end if
    end do

    do i = 1, size(enlist, 2)
      ele_owner(i) = minval(node_owner(enlist(:, i)))
    end do

    deallocate(node_owner)

  end subroutine element_ownership

  subroutine node_ownership(nnodes, enlist, halo, node_owner)
    integer, intent(in) :: nnodes
    ! loc x nelements
    integer, dimension(:, :), intent(in) :: enlist
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

end module libsupermesh_halo_ownership
