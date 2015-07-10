#include "fdebug.h"

module libsupermesh_construction

  use libsupermesh_fldebug
  use libsupermesh_fields_dummy
  use libsupermesh_tri_intersection_module
  use libsupermesh_tet_intersection_module
  
  implicit none

  public :: intersect_elements, intersect_elements_old
  
  private
 
  contains
  
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
    type(mesh_type) :: intersection_mesh

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
      call libsupermesh_cintersector_get_output(nonods, totele, ndimA, locA, intersection%val, intersection_mesh%ndglno)
    end if

    call deallocate(intersection_mesh)
    
  end function intersect_elements_old
  
end module libsupermesh_construction
