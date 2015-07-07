#include "fdebug.h"

module libsupermesh_construction

  use libsupermesh_fldebug
  use libsupermesh_fields_dummy
  use libsupermesh_global_parameters, only : real_8
  use libsupermesh_tri_intersection_module
  use libsupermesh_tet_intersection_module
  
  implicit none
  
  interface libsupermesh_cintersector_set_input  
    subroutine libsupermesh_cintersector_set_input(nodes_A, nodes_B, ndim, loc)
      use libsupermesh_global_parameters, only : real_8
      implicit none
      real(kind = real_8), dimension(ndim, loc), intent(in) :: nodes_A, nodes_B
      integer, intent(in) :: ndim, loc
    end subroutine libsupermesh_cintersector_set_input
  end interface libsupermesh_cintersector_set_input
  
  interface 
    subroutine libsupermesh_cintersector_drive
    end subroutine libsupermesh_cintersector_drive
  end interface

  interface
    subroutine libsupermesh_cintersector_query(nonods, totele)
      implicit none
      integer, intent(out) :: nonods, totele
    end subroutine libsupermesh_cintersector_query
  end interface

  interface libsupermesh_cintersector_get_output  
    subroutine libsupermesh_cintersector_get_output(nonods, totele, ndim, loc, nodes, enlist)
      use libsupermesh_global_parameters, only : real_8
      implicit none
      integer, intent(in) :: nonods, totele, ndim, loc
      real(kind = real_8), dimension(nonods * ndim), intent(out) :: nodes
      integer, dimension(totele * loc), intent(out) :: enlist
    end subroutine libsupermesh_cintersector_get_output
  end interface libsupermesh_cintersector_get_output
  
  interface intersector_set_dimension
    subroutine libsupermesh_cintersector_set_dimension(ndim)
      implicit none
      integer, intent(in) :: ndim
    end subroutine libsupermesh_cintersector_set_dimension
  end interface intersector_set_dimension
  
  interface 
    subroutine libsupermesh_cintersector_set_exactness(exact)
      implicit none
      integer, intent(in) :: exact
    end subroutine libsupermesh_cintersector_set_exactness
  end interface
  
  ! I hope this is big enough ...
  real, dimension(1024) :: nodes_tmp
  logical :: intersector_exactness = .false.
  type(mesh_type), save :: intersection_mesh
  logical, save :: mesh_allocated = .false.

  public :: intersect_elements, intersect_elements_old, &
          &    intersector_set_dimension, &
          &    intersector_set_exactness, finalise_libsupermesh
  public :: intersector_exactness
  
  private
 
  contains
  
  subroutine finalise_libsupermesh
    if (mesh_allocated) then
      call deallocate(intersection_mesh)
      mesh_allocated = .false.
    end if
  end subroutine finalise_libsupermesh
  
  subroutine intersect_elements(positions_A_val, &
        posB, n_C, simplexC_real)

    real, dimension(:, :), intent(in) :: positions_A_val
    real, dimension(:, :), intent(in) :: posB

    real, dimension(size(positions_A_val, 1), size(positions_A_val, 1) + 1, simplex_size(size(positions_A_val, 1))), intent(inout) :: simplexC_real

    integer, intent(out) :: n_C

    integer :: ndimA
    type(tri_type) :: triA, triB
    type(tet_type) :: tetA, tetB

    ewrite(1, *) "In intersect_elements"

    ndimA = size(positions_A_val, 1)
    select case(ndimA)
    case(1)
      call intersect_intervals(positions_A_val(1, :), posB(1, :), simplexC_real(1, :, 1), n_C)
    case(2)
      triA%v = positions_A_val
      triB%v = posB

      call intersect_tris(triA%v, triB%v, simplexC_real, n_C)
    case(3)
      tetA%v = positions_A_val
      tetB%v = posB

      call intersect_tets(tetA%v, tetB%v, simplexC_real, n_C)
    case default
      FLAbort("intersect_elements: Unsupported dimension.")
    end select

  end subroutine intersect_elements

  pure function simplex_size(dim) result(size)
    integer, intent(in) :: dim

    integer :: size

    select case(dim)
      case(1)
        size = 1
      case(2)
        size = tri_buf_size
      case(3)
        size = tet_buf_size
      case default
        size = 0
    end select

  end function simplex_size
  
  subroutine intersect_intervals(positions_A_val, posB, intC, n_C)
    real, dimension(2), intent(in) :: positions_A_val
    real, dimension(2), intent(in) :: posB
    real, dimension(2), intent(out) :: intC
    integer, intent(out) :: n_C

    real :: min_a, max_a, min_b, max_b

    min_a = minval(positions_A_val)
    max_a = maxval(positions_A_val)
    min_b = minval(posB)
    max_b = maxval(posB)

    if(max_b <= min_a .OR. min_b >= max_a) then
      n_C = 0
    else
      intC(1) = max(min_a, min_b)
      intC(2) = min(max_a, max_b)
      n_C = 1
    end if

  end subroutine intersect_intervals

  function intersect_elements_old(positions_A_val, &
        posB, locA, ndimA) result(intersection)

    real, intent(in), dimension(ndimA, locA) :: positions_A_val
    integer, intent(in) :: locA, ndimA
    type(vector_field) :: intersection
    real, dimension(:, :), intent(in) :: posB
  
    integer :: i, nonods, totele

    call libsupermesh_cintersector_set_input(positions_A_val, posB, ndimA, locA)
    call libsupermesh_cintersector_drive
    call libsupermesh_cintersector_query(nonods, totele)
    call allocate(intersection_mesh, ndimA, nonods, totele, ndimA + 1)
    intersection_mesh%ndglno = (/ (i, i=1,totele) /)
    intersection_mesh%continuity = -1
    call allocate(intersection, ndimA, intersection_mesh)
    if (nonods > 0) then
#ifdef DDEBUG
      intersection_mesh%ndglno = -1
#endif
      call libsupermesh_cintersector_get_output(nonods, totele, ndimA, locA, nodes_tmp, intersection_mesh%ndglno)

      do i = 1, ndimA
        intersection%val(i,:) = nodes_tmp((i - 1) * nonods + 1:i * nonods)
      end do
    end if

    call deallocate(intersection_mesh)
    
  end function intersect_elements_old
  
  subroutine intersector_set_exactness(exactness)
    logical, intent(in) :: exactness
    integer :: exact

    if (exactness) then
      exact = 1
    else
      exact = 0
    end if
    intersector_exactness = exactness

    call libsupermesh_cintersector_set_exactness(exact)
  end subroutine intersector_set_exactness
  
end module libsupermesh_construction
