#include "libsupermesh_debug.h"

module libsupermesh_intersections

  implicit none
  
  private
  
  public :: intersections, deallocate, intersections_to_csr_sparsity
  
  type intersections
    integer, dimension(:), pointer :: v
    integer :: n
  end type intersections

  interface deallocate
    module procedure deallocate_intersections, deallocate_intersections_rank_1
  end interface deallocate
  
contains

  pure subroutine deallocate_intersections(ints)
    type(intersections), intent(inout) :: ints

    deallocate(ints%v)

  end subroutine deallocate_intersections

  pure subroutine deallocate_intersections_rank_1(ints)
    type(intersections), dimension(:), intent(inout) :: ints

    integer :: i

    do i = 1, size(ints)
      deallocate(ints(i)%v)
    end do

  end subroutine deallocate_intersections_rank_1

  pure subroutine intersections_to_csr_sparsity(ints, indices, indptr)
    type(intersections), dimension(:), intent(in) :: ints
    ! Compressed Sparse Row (CSR) sparsity pattern, as described in:
    !   http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.sparse.csr_matrix.html
    integer, dimension(:), allocatable, intent(out) :: indices
    integer, dimension(:), intent(out) :: indptr

    integer :: i, j, n

    n = 0
    do i = 1, size(ints)
      n = n + ints(i)%n
    end do

    allocate(indices(n))
    indptr(1) = 1
    do i = 1, size(ints)
      indptr(i + 1) = indptr(i) + ints(i)%n
      do j = 1, ints(i)%n
        indices(indptr(i) + j - 1) = ints(i)%v(j)
      end do
    end do

  end subroutine intersections_to_csr_sparsity
  
end module libsupermesh_intersections
