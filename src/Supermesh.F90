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
  
  function libsupermesh_intersect_elements(positions_A_val, elementCountA, verticesA, quadDimA, ele_A, &
        posB, shape, locA, ndimA, nnodesA, fieldMeshShapeLocA, fieldTypeA, &
        positions_a_MeshNdglno) result(intersection)
    real, intent(in), dimension(ndimA, nnodesA) :: positions_A_val
    integer, intent(in), dimension(elementCountA * fieldMeshShapeLocA) :: positions_a_MeshNdglno
    integer, intent(in) :: ele_A, locA, ndimA, nnodesA, fieldMeshShapeLocA, fieldTypeA
    type(vector_field) :: intersection
    type(mesh_type) :: intersection_mesh
    type(element_type), intent(in) :: shape
    real, dimension(:, :), intent(in) :: posB
  
    integer :: dim, elementCountA, verticesA, quadDimA
    integer :: nonods, totele
    integer :: i
    
    type(vector_field), target :: positions_A
    type(mesh_type) :: mesh_lib
    type(quadrature_type) :: quad_lib
    type(element_type) :: shape_lib

    ewrite(1, *) "In libsupermesh_intersect_elements"

#ifdef DDEBUG
    select case(ndimA)
      case(2)
        assert(shape%loc == 3)
      case(3)
        assert(shape%loc == 4)
    end select
#endif

    quad_lib = make_quadrature(vertices = verticesA, dim = quadDimA, ngi = 1, degree = 2)
    shape_lib = make_element_shape(vertices = fieldMeshShapeLocA, dim = ndimA, degree = 1, quad = quad_lib)
    call allocate(mesh_lib, nnodesA, elementCountA, shape_lib)
    
    mesh_lib%ndglno = positions_a_MeshNdglno
    call allocate(positions_A, ndimA, mesh_lib)
    positions_A%val = positions_A_val 
    positions_A%dim = ndimA

    call libsupermesh_cintersector_set_input(ele_val(positions_A, ele_A), posB, ndimA, locA)
    call libsupermesh_cintersector_drive
    call libsupermesh_cintersector_query(nonods, totele)
    ! IAKOVOS REMOVE COMMENT
!    write(*,*) "libsupermesh_intersect_elements: nonods:",nonods,", totele:",totele,"."
    call allocate(intersection_mesh, nonods, totele, shape, "IntersectionMesh")
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

    call deallocate(quad_lib)
    call deallocate(shape_lib)
    call deallocate(mesh_lib)
    
    call deallocate(intersection_mesh)
    call deallocate(positions_A)

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
