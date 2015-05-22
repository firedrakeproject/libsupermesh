#include "fdebug.h"

module libsupermesh_construction
  use libsupermesh_fields_data_types
  use libsupermesh_fields_allocates
  use libsupermesh_fields_base
  use libsupermesh_sparse_tools
  use libsupermesh_futils
!  use metric_tools			! IAKOVOS commented out
  use libsupermesh_linked_lists
!  use unify_meshes_module		! IAKOVOS commented out
  use libsupermesh_transform_elements
  use libsupermesh_global_parameters, only : real_4, real_8
!  use tetrahedron_intersection_module	! IAKOVOS commented out
  use libsupermesh_tet_intersection_module
  use libsupermesh_shape_functions, only: make_element_shape
  implicit none
  
  interface libsupermesh_cintersector_set_input
    module procedure libsupermesh_intersector_set_input_sp
  
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
    module procedure libsupermesh_intersector_get_output_sp
  
    subroutine libsupermesh_cintersector_get_output(nonods, totele, ndim, loc, nodes, enlist)
      use libsupermesh_global_parameters, only : real_8
      implicit none
      integer, intent(in) :: nonods, totele, ndim, loc
      real(kind = real_8), dimension(nonods * ndim), intent(out) :: nodes
      integer, dimension(totele * loc), intent(out) :: enlist
    end subroutine libsupermesh_cintersector_get_output
  end interface libsupermesh_cintersector_get_output
  
  interface libsupermesh_intersector_set_dimension
    subroutine libsupermesh_cintersector_set_dimension(ndim)
      implicit none
      integer, intent(in) :: ndim
    end subroutine libsupermesh_cintersector_set_dimension
  end interface libsupermesh_intersector_set_dimension
  
  interface 
    subroutine libsupermesh_cintersector_set_exactness(exact)
      implicit none
      integer, intent(in) :: exact
    end subroutine libsupermesh_cintersector_set_exactness
  end interface
  
  ! I hope this is big enough ...
  real, dimension(1024) :: nodes_tmp
  logical :: libsupermesh_intersector_exactness = .false.

  public :: libsupermesh_intersect_elements, libsupermesh_intersector_set_dimension, &
                  libsupermesh_intersector_set_exactness
  public :: libsupermesh_intersector_exactness

  private
 
  contains
  
  subroutine libsupermesh_intersector_set_input_sp(nodes_A, nodes_B, ndim, loc)
    real(kind = real_4), dimension(ndim, loc), intent(in) :: nodes_A
    real(kind = real_4), dimension(ndim, loc), intent(in) :: nodes_B
    integer, intent(in) :: ndim
    integer, intent(in) :: loc
    
    call libsupermesh_cintersector_set_input(real(nodes_A, kind = real_8), real(nodes_B, kind = real_8), ndim, loc)
  
  end subroutine libsupermesh_intersector_set_input_sp
  
  subroutine libsupermesh_intersector_get_output_sp(nonods, totele, ndim, loc, nodes, enlist)
    integer, intent(in) :: nonods
    integer, intent(in) :: totele
    integer, intent(in) :: ndim
    integer, intent(in) :: loc
    real(kind = real_4), dimension(nonods * ndim), intent(out) :: nodes
    integer, dimension(totele * loc), intent(out) :: enlist
    
    real(kind = real_8), dimension(size(nodes)) :: lnodes
    
    call libsupermesh_cintersector_get_output(nonods, totele, ndim, loc, lnodes, enlist)
    nodes = lnodes
    
  end subroutine libsupermesh_intersector_get_output_sp
  
  function libsupermesh_intersect_elements(positions_A_val, &
        posB, locA, ndimA, nnodesA, &
        quadVertices, quadDim, quadNgi, quadDegree, &
        shapeLoc, shapeDim, shapeDegree) result(intersection)
    real, intent(in), dimension(ndimA, locA) :: positions_A_val
    integer, intent(in) :: locA, ndimA, nnodesA
    type(vector_field) :: intersection
    type(mesh_type) :: intersection_mesh
    real, dimension(:, :), intent(in) :: posB
  
    integer :: nonods, totele, i
    integer, intent(in) :: quadVertices, quadDim, quadNgi, quadDegree, shapeLoc, shapeDim, shapeDegree

    ewrite(1, *) "In libsupermesh_intersect_elements"
#ifdef DDEBUG
    select case(ndimA)
      case(2)
        assert(shapeLoc == 3)
      case(3)
        assert(shapeLoc == 4)
    end select
#endif

    call libsupermesh_cintersector_set_input(positions_A_val, posB, ndimA, locA)
    call libsupermesh_cintersector_drive
    call libsupermesh_cintersector_query(nonods, totele)
    call allocate(intersection_mesh, nonods, totele, shapeLoc, shapeDim, shapeDegree, quadVertices, quadDim, quadNgi, quadDegree, name="IntersectionMesh")
!    intersection_mesh%ndglno = (/ (i, i=1,totele) /)
    intersection_mesh%continuity = -1
    call allocate(intersection, ndimA, intersection_mesh, "IntersectionCoordinates")
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

  end function libsupermesh_intersect_elements
  
  subroutine libsupermesh_intersector_set_exactness(exactness)
    logical, intent(in) :: exactness
    integer :: exact

    if (exactness) then
      exact = 1
    else
      exact = 0
    end if
    libsupermesh_intersector_exactness = exactness

    call libsupermesh_cintersector_set_exactness(exact)
  end subroutine libsupermesh_intersector_set_exactness
  
end module libsupermesh_construction
