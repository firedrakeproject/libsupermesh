#include "fdebug.h"

module libsupermesh_parallel_supermesh

  use libsupermesh_integer_hash_table
  use libsupermesh_fldebug
  use libsupermesh_read_halos

  implicit none

#include "mpif.h"

contains

  subroutine parallel_supermesh(positions_a, enlist_a, un_a, halo_a, positions_b, enlist_b, halo_b, donor_ele_data, unpack_donor_ele_data, intersection_calculation)
    ! dim x nnodes_a
    real, dimension(:, :), intent(in)    :: positions_a
    ! loc_a x nelements_a
    integer, dimension(:, :), intent(in) :: enlist_a
    ! nnodes_a
    integer, dimension(:), intent(in)    :: un_a
    type(halo_type), intent(in)          :: halo_a
    ! dim x nnodes_b
    real, dimension(:, :), intent(in)    :: positions_b
    ! loc_b x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    type(halo_type), intent(in)          :: halo_b
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

    integer                            :: parallel_ele_A, parallel_ele_B, rank, ierr
    integer, dimension(:), allocatable :: ele_owner_a, ele_owner_b
    real, dimension(:,:), allocatable  :: bbox_a, bbox_b

!     call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
!     assert(ierr == MPI_SUCCESS)

    ! 1. Communicate bounding box data for donor and target
    ! 
!     allocate(ele_owner_a(size(enlist_a(:,1))))
!     call element_ownership(node_count(positions_a), enlist_a, halo_a, ele_owner_a)
!     parallel_ele_A = count(ele_owner_a == rank)
! 
!     allocate(ele_owner_b(ele_count(positions_b)))
!     call element_ownership(node_count(positions_b), enlist_b, halo_b, ele_owner_b)
!     parallel_ele_B = count(ele_owner_b == rank)
! 
!     allocate(bbox_a(positions_a,2))
!     bbox_a = partition_bbox(positions_a, ele_owner_a, rank)
!     allocate(bbox_b(positions_b,2))
!     bbox_b = partition_bbox(positions_b, ele_owner_b, rank)
    ! 2. Use bounding box data to cull donor mesh

    ! 3. Pack non-culled mesh data for communication, calling user specified data functions

    ! 4. Communicate donor mesh and mesh data

    ! 5. Supermesh and call user specified element unpack and calculation functions
    
  end subroutine parallel_supermesh

end module libsupermesh_parallel_supermesh
