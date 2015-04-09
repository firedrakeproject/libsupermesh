#include "fdebug.h"

module libsupermesh_construction
  use libsupermesh_fields_data_types, mesh_faces_lib => mesh_faces, &
                    mesh_subdomain_mesh_lib => mesh_subdomain_mesh, &
                    scalar_field_lib => scalar_field, &
                    vector_field_lib => vector_field, &
                    tensor_field_lib => tensor_field, &
  scalar_boundary_condition_lib => scalar_boundary_condition, &
  vector_boundary_condition_lib => vector_boundary_condition, &
  scalar_boundary_conditions_ptr_lib => scalar_boundary_conditions_ptr, &
  vector_boundary_conditions_ptr => vector_boundary_conditions_ptr
  use libsupermesh_fields_allocates
  use libsupermesh_fields_base
  use libsupermesh_sparse_tools, wrap_lib => wrap
  use libsupermesh_futils, real_format_lib => real_format, &
        real_format_len_lib => real_format_len
!  use metric_tools			! IAKOVOS commented out
  use libsupermesh_linked_lists
!  use unify_meshes_module		! IAKOVOS commented out
  use libsupermesh_transform_elements
  use libsupermesh_global_parameters, only : real_4, real_8
!  use tetrahedron_intersection_module	! IAKOVOS commented out
  use libsupermesh_tet_intersection_module
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
  
  function libsupermesh_intersect_elements(positions_A, ele_A, posB, shape, locA, ndimA, nnodesA, fieldMeshShapeLocA, fieldTypeA, positions_a_MeshNdglno) result(intersection)
    real, intent(in), dimension(ndimA, fieldMeshShapeLocA) :: positions_A
    integer, intent(in), dimension(nnodesA * ndimA) :: positions_a_MeshNdglno
    integer, intent(in) :: ele_A, locA, ndimA, nnodesA, fieldMeshShapeLocA, fieldTypeA
    type(vector_field_lib) :: intersection
    type(mesh_type) :: intersection_mesh
    type(element_type), intent(in) :: shape
    real, dimension(:, :), intent(in) :: posB
  
    integer :: dim
    integer :: nonods, totele
    integer :: i

!    dim = positions_A%dim
#ifdef DDEBUG
    select case(ndimA)
      case(2)
        assert(shape%loc == 3)
      case(3)
        assert(shape%loc == 4)
    end select
#endif

    call libsupermesh_cintersector_set_input(ele_val_v(positions_A, ele_A, ndimA, nnodesA, fieldMeshShapeLocA, fieldTypeA, positions_a_MeshNdglno), posB, ndimA, locA)
    call libsupermesh_cintersector_drive
    call libsupermesh_cintersector_query(nonods, totele)
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
