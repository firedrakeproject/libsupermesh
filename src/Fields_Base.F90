!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineeringp
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA
#include "fdebug.h"

module libsupermesh_fields_base
  !!< This module contains abstracted field types which carry shape and
  !!< connectivity with them.
  use libsupermesh_shape_functions, only: element_type
!  use tensors					! IAKOVOS commented out
  use libsupermesh_fields_data_types
  use libsupermesh_reference_counting
  use libsupermesh_global_parameters, only: FIELD_NAME_LEN, current_debug_level, current_time
  use libsupermesh_elements
  use libsupermesh_element_numbering
!  use embed_python				! IAKOVOS commented out
  use libsupermesh_sparse_tools
  use libsupermesh_vector_tools, only: solve, invert, norm2, cross_product
  implicit none
  
  interface ele_nodes
     module procedure ele_nodes_scalar, ele_nodes_vector, ele_nodes_tensor,&
          & ele_nodes_mesh
  end interface
  
  interface face_local_nodes
!     module procedure face_local_nodes_mesh,face_local_nodes_scalar,&
!          & face_local_nodes_vector, face_local_nodes_tensor
     module procedure face_local_nodes_mesh
  end interface
  
  interface face_global_nodes
     module procedure face_global_nodes_mesh, face_global_nodes_vector,&
          & face_global_nodes_scalar, face_global_nodes_tensor
  end interface
  
  interface ele_neigh
!     module procedure ele_neigh_mesh, ele_neigh_scalar, ele_neigh_vector, &
!          & ele_neigh_tensor, ele_neigh_i_mesh, ele_neigh_i_scalar, &
!          & ele_neigh_i_vector, ele_neigh_i_tensor
     module procedure ele_neigh_mesh, ele_neigh_vector
  end interface
  
  interface ele_faces
!     module procedure ele_faces_mesh, ele_faces_vector, ele_faces_scalar, &
!          & ele_faces_tensor
    module procedure ele_faces_mesh, ele_faces_vector
  end interface
  
  interface ele_face
!     module procedure ele_face_mesh, ele_face_scalar, ele_face_vector,&
!          & ele_face_tensor
     module procedure ele_face_mesh, ele_face_vector
  end interface
  
  interface ele_shape
!     module procedure ele_shape_scalar, ele_shape_vector, ele_shape_tensor,&
!          & ele_shape_mesh 
     module procedure ele_shape_vector, ele_shape_mesh
  end interface
  
  interface ele_ngi
!     module procedure ele_ngi_scalar, ele_ngi_vector, ele_ngi_tensor,&
!          & ele_ngi_mesh
     module procedure ele_ngi_vector
  end interface
  
  interface face_loc
     module procedure face_loc_scalar, face_loc_vector, face_loc_tensor,&
          & face_loc_mesh
  end interface
  
  interface face_ele
!     module procedure face_ele_mesh, face_ele_scalar, face_ele_vector,&
!          & face_ele_tensor
      module procedure face_ele_mesh, face_ele_vector
  end interface
  
  interface face_ngi
!     module procedure face_ngi_scalar, face_ngi_vector, face_ngi_tensor,&
!          & face_ngi_mesh 
     module procedure face_ngi_vector
  end interface
  
  interface face_shape
     module procedure face_shape_scalar, face_shape_vector,&
          & face_shape_tensor, face_shape_mesh
  end interface
  
  interface mesh_dim
     module procedure mesh_dim_mesh, mesh_dim_scalar, mesh_dim_vector,&
          & mesh_dim_tensor
  end interface
  
  interface extract_scalar_field ! extract_scalar_field is already used in State.F90
!     module procedure extract_scalar_field_from_vector_field, extract_scalar_field_from_tensor_field
     module procedure extract_scalar_field_from_vector_field
  end interface
  
  interface halo_count
    module procedure halo_count_mesh, halo_count_scalar, halo_count_vector, &
      & halo_count_tensor
  end interface halo_count
  
  interface element_halo_count
    module procedure element_halo_count_mesh, element_halo_count_scalar, &
      & element_halo_count_vector, element_halo_count_tensor
  end interface element_halo_count
  
  interface operator(==)
     module procedure mesh_equal
  end interface
  
  interface surface_element_count
!    module procedure surface_element_count_scalar, surface_element_count_vector, &
!          & surface_element_count_tensor, surface_element_count_mesh
     module procedure surface_element_count_vector, surface_element_count_mesh
  end interface
  
  interface set_from_python_function
     module procedure set_values_from_python_scalar, set_values_from_python_scalar_pos, &
       set_values_from_python_vector, set_values_from_python_vector_pos, &
       set_values_from_python_vector_field
  end interface

  interface ele_loc
     module procedure ele_loc_scalar, ele_loc_vector, ele_loc_tensor,&
          & ele_loc_mesh
  end interface
  
  interface face_vertices
!     module procedure face_vertices_scalar, face_vertices_vector, face_vertices_tensor,&
!          & face_vertices_mesh, face_vertices_shape
    module procedure face_vertices_mesh, face_vertices_shape
  end interface
  
  interface ele_val
     module procedure ele_val_scalar, ele_val_vector, ele_val_vector_dim,&
          & ele_val_tensor, ele_val_tensor_dim_dim
  end interface
  
  interface face_val
     module procedure face_val_scalar, face_val_vector, face_val_tensor,&
          & face_val_vector_dim, face_val_tensor_dim_dim
  end interface
  
  interface has_discontinuous_internal_boundaries
     module procedure mesh_has_discontinuous_internal_boundaries
  end interface has_discontinuous_internal_boundaries

  interface element_count
     module procedure element_count_scalar, element_count_vector,&
          & element_count_tensor, element_count_mesh
  end interface
  
  interface ele_count
     module procedure element_count_scalar, element_count_vector,&
          & element_count_tensor, element_count_mesh
  end interface
  
  interface face_count
    module procedure face_count_scalar, face_count_vector, &
          & face_count_tensor, face_count_mesh
  end interface
  
  interface node_count
     module procedure node_count_scalar, node_count_vector,&
          & node_count_tensor, node_count_mesh
  end interface
  
  interface local_coords
!    module procedure local_coords_interpolation, &
!          local_coords_interpolation_all, local_coords_mesh, &
!          local_coords_scalar, local_coords_vector, local_coords_tensor
    module procedure local_coords_interpolation, &
           local_coords_vector, local_coords_mesh
  end interface
  
  interface local_face_number
!     module procedure local_face_number_mesh, local_face_number_scalar, &
!          & local_face_number_vector, local_face_number_tensor
    module procedure local_face_number_mesh
  end interface local_face_number
  
  interface node_val
     module procedure node_val_scalar, node_val_vector, node_val_tensor, &
          & node_val_scalar_v, node_val_vector_v, node_val_vector_dim_v,&
          & node_val_tensor_v, node_val_vector_dim, node_val_tensor_dim_dim, &
          & node_val_tensor_dim_dim_v
  end interface
  
  interface has_faces
     module procedure has_faces_mesh
  end interface
  
  interface tetvol
    module procedure tetvol_old
  end interface tetvol
  
contains

  pure function mesh_dim_mesh(mesh) result (mesh_dim)
    ! Return the dimensionality of the mesh.
    integer :: mesh_dim
    type(mesh_type), intent(in) :: mesh
    
    mesh_dim=mesh%shape%dim

  end function mesh_dim_mesh
  
  pure function mesh_dim_scalar(field) result (mesh_dim)
    ! Return the dimensionality of the field.
    integer :: mesh_dim
    type(scalar_field), intent(in) :: field
    
    mesh_dim=field%mesh%shape%dim

  end function mesh_dim_scalar

  pure function mesh_dim_vector(field) result (mesh_dim)
    ! Return the dimensionality of the field.
    integer :: mesh_dim
    type(vector_field), intent(in) :: field
    
    mesh_dim=field%mesh%shape%dim

  end function mesh_dim_vector

  pure function mesh_dim_tensor(field) result (mesh_dim)
    ! Return the dimensionality of the field.
    integer :: mesh_dim
    type(tensor_field), intent(in) :: field
    
    mesh_dim=field%mesh%shape%dim

  end function mesh_dim_tensor
  
  function mesh_connectivity(mesh) result (ndglno)
    !!< Assuming that the input mesh is at least C0, return the connectivity
    !!< of the mesh.
    type(mesh_type), intent(in) :: mesh
    integer, dimension(mesh%elements*mesh%shape%numbering%vertices) ::&
         & ndglno

    integer, dimension(mesh%shape%numbering%vertices) :: vertices
    integer :: i, nodes

    integer, dimension(:), allocatable :: map
    integer, dimension(:), pointer :: e_nodes

    vertices=local_vertices(mesh%shape%numbering)

    do i=1,mesh%elements
       e_nodes => ele_nodes(mesh, i)
       ndglno((i-1)*size(vertices)+1:i*size(vertices)) = e_nodes(vertices)
    end do

    allocate(map(maxval(ndglno)))

    map=0
    nodes=0

    do i=1,size(ndglno)
       if (map(ndglno(i))==0) then
          nodes=nodes+1
          map(ndglno(i))=nodes
       end if
       
       ndglno(i)=map(ndglno(i))
    end do
    
  end function mesh_connectivity
  
  pure function mesh_equal(mesh1, mesh2)
    !!< Test for the equality of two meshes. This is not totally safe since
    !!< we do not compare all of ndglno. We assume that two meshes are equal
    !!< if they have the same element, continuity, node count and element
    !!< count. This should be sufficient in almost all circumstances.
    logical :: mesh_equal
    type(mesh_type), intent(in) :: mesh1, mesh2

    if (associated(mesh1%refcount) .and. associated(mesh2%refcount)) then
       mesh_equal= (mesh1%refcount%id==mesh2%refcount%id)
    else
       mesh_equal= .false.
    end if

  end function mesh_equal
  
  pure function halo_count_mesh(mesh) result(count)
    type(mesh_type), intent(in) :: mesh
    
    integer :: count
    
    if(.not. associated(mesh%halos)) then
      count = 0
    else
      count = size(mesh%halos)
    end if
  
  end function halo_count_mesh
  
  pure function halo_count_scalar(s_field) result(count)
    type(scalar_field), intent(in) :: s_field
    
    integer :: count
    
    count = halo_count(s_field%mesh)
  
  end function halo_count_scalar
  
  pure function halo_count_vector(v_field) result(count)
    type(vector_field), intent(in) :: v_field
    
    integer :: count
    
    count = halo_count(v_field%mesh)
  
  end function halo_count_vector
  
  pure function halo_count_tensor(t_field) result(count)
    type(tensor_field), intent(in) :: t_field
    
    integer :: count
    
    count = halo_count(t_field%mesh)
  
  end function halo_count_tensor

  pure function element_halo_count_mesh(mesh) result(count)
    type(mesh_type), intent(in) :: mesh
    
    integer :: count
    
    if(.not. associated(mesh%element_halos)) then
      count = 0
    else
      count = size(mesh%element_halos)
    end if
  
  end function element_halo_count_mesh
  
  pure function element_halo_count_scalar(s_field) result(count)
    type(scalar_field), intent(in) :: s_field
    
    integer :: count
    
    count = element_halo_count(s_field%mesh)
  
  end function element_halo_count_scalar
  
  pure function element_halo_count_vector(v_field) result(count)
    type(vector_field), intent(in) :: v_field
    
    integer :: count
    
    count = element_halo_count(v_field%mesh)
  
  end function element_halo_count_vector
  
  pure function element_halo_count_tensor(t_field) result(count)
    type(tensor_field), intent(in) :: t_field
    
    integer :: count
    
    count = element_halo_count(t_field%mesh)
  
  end function element_halo_count_tensor

  function local_coords_interpolation(position_field, ele, position) result(local_coords)
    !!< Given a position field, this returns the local coordinates of
    !!< position with respect to element "ele".
    !!<
    !!< This assumes the position field is linear. For higher order
    !!< only the coordinates of the vertices are considered
    type(vector_field), intent(in) :: position_field
    integer, intent(in) :: ele
    real, dimension(:), intent(in) :: position
    real, dimension(size(position) + 1) :: local_coords
    real, dimension(mesh_dim(position_field) + 1, size(position) + 1) :: matrix
    real, dimension(mesh_dim(position_field), size(position) + 1) :: tmp_matrix
    integer, dimension(position_field%mesh%shape%numbering%vertices):: vertices
    integer, dimension(:), pointer:: nodes
    integer :: dim

    dim = size(position)

    assert(dim == mesh_dim(position_field))
    assert(position_field%mesh%shape%numbering%family==FAMILY_SIMPLEX)
    assert(position_field%mesh%shape%numbering%type==ELEMENT_LAGRANGIAN)

    local_coords(1:dim) = position
    local_coords(dim+1) = 1.0

    if (position_field%mesh%shape%degree==1) then
      tmp_matrix = ele_val(position_field, ele)
    else
      nodes => ele_nodes(position_field, ele)
      vertices=local_vertices(position_field%mesh%shape%numbering)
      tmp_matrix = node_val(position_field, nodes(vertices) )
    end if
    
    matrix(1:dim, :) = tmp_matrix
    matrix(dim+1, :) = 1.0

    call solve(matrix, local_coords)
    
  end function local_coords_interpolation
  
  function local_coords_vector(field, ele, node, stat) result(local_coord)
    !!< returns the local node number within a given element ele of the global
    !!< node number node
    type(vector_field), intent(in) :: field
    integer, intent(in) :: ele, node
    integer, intent(inout), optional :: stat
    integer :: local_coord
    
    local_coord = local_coords(field%mesh, ele, node, stat)
    
  end function local_coords_vector
  
  function local_coords_mesh(mesh, ele, node, stat) result(local_coord)
    !!< returns the local node number within a given element ele of the global
    !!< node number node
    type(mesh_type), intent(in) :: mesh
    integer, intent(in) :: ele, node
    integer, intent(inout), optional :: stat
    integer :: local_coord
    integer, dimension(:), pointer :: nodes
    integer :: i

    if(present(stat)) stat = 0

    local_coord = 0

    nodes => ele_nodes(mesh, ele)
    do i=1,size(nodes)
      if (nodes(i) == node) then
        local_coord = i
        return
      end if
    end do
    
    if(present(stat)) then
      stat = 1
    else
      FLAbort("Failed to find local coordinate.")
    end if
    
  end function local_coords_mesh
  
  function ele_faces_mesh(mesh, ele_number) result (ele_faces)
    !!< Return a pointer to a vector containing the face numbers of the
    !!< faces adjacent to ele_number.
    integer, dimension(:), pointer :: ele_faces
    type(mesh_type),intent(in) :: mesh
    integer, intent(in) :: ele_number

    ele_faces =>row_ival_ptr(mesh%faces%face_list, ele_number)
  
  end function ele_faces_mesh
  
  function ele_neigh_mesh(mesh, ele_number) result (ele_neigh)
    !!< Return a pointer to a vector containing the element numbers of the
    !!< elements adjacent to ele_number.
    integer, dimension(:), pointer :: ele_neigh
    type(mesh_type),intent(in) :: mesh
    integer, intent(in) :: ele_number

    ele_neigh=>row_m_ptr(mesh%faces%face_list, ele_number)
  
  end function ele_neigh_mesh
  
  function ele_neigh_vector(field, ele_number) result (ele_neigh)
    !!< Return a pointer to a vector containing the element numbers of the
    !!< elements adjacent to ele_number.
    integer, dimension(:), pointer :: ele_neigh
    type(vector_field),intent(in) :: field
    integer, intent(in) :: ele_number

    ele_neigh=>ele_neigh_mesh(field%mesh, ele_number)
  
  end function ele_neigh_vector
  
  function ele_faces_vector(field, ele_number) result (ele_faces)
    !!< Return a pointer to a vector containing the face numbers of the
    !!< faces adjacent to ele_number.
    integer, dimension(:), pointer :: ele_faces
    type(vector_field),intent(in) :: field
    integer, intent(in) :: ele_number

    ele_faces=>ele_faces_mesh(field%mesh, ele_number)
  
  end function ele_faces_vector
  
  pure function mesh_has_discontinuous_internal_boundaries(mesh)
    !!< Return whether the mesh has discontinuous boundaries
    !!< These are internal boundaries where the surface id are 
    !!< allowed to be discontinuous (the pair of adjacent 
    !!< internal facets can have two different ids).
    logical :: mesh_has_discontinuous_internal_boundaries
    type(mesh_type), intent(in) :: mesh
    
    if (associated(mesh%faces)) then
      mesh_has_discontinuous_internal_boundaries = mesh%faces%has_discontinuous_internal_boundaries
    else
      mesh_has_discontinuous_internal_boundaries = .false.
    end if
  
  end function mesh_has_discontinuous_internal_boundaries

  pure function ele_ngi_vector(field, ele_number) result (ele_ngi)
    ! Return the number of nodes of element ele_number.
    integer :: ele_ngi
    type(vector_field),intent(in) :: field
    integer, intent(in) :: ele_number
    
    ele_ngi=field%mesh%shape%ngi
    
  end function ele_ngi_vector

  pure function face_loc_mesh(mesh, face_number) result (face_loc)
    ! Return the number of nodes of face face_number.
    integer :: face_loc
    type(mesh_type),intent(in) :: mesh
    integer, intent(in) :: face_number
    
    face_loc=mesh%faces%shape%loc

  end function face_loc_mesh
  
  pure function face_loc_scalar(field, face_number) result (face_loc)
    ! Return the number of nodes of face face_number.
    integer :: face_loc
    type(scalar_field),intent(in) :: field
    integer, intent(in) :: face_number
    
    face_loc=field%mesh%faces%shape%loc
    
  end function face_loc_scalar
  
  pure function face_loc_vector(field, face_number) result (face_loc)
    ! Return the number of nodes of face face_number.
    integer :: face_loc
    type(vector_field),intent(in) :: field
    integer, intent(in) :: face_number
    
    face_loc=field%mesh%faces%shape%loc
    
  end function face_loc_vector
  
  pure function face_loc_tensor(field, face_number) result (face_loc)
    ! Return the number of nodes of face face_number.
    integer :: face_loc
    type(tensor_field),intent(in) :: field
    integer, intent(in) :: face_number
    
    face_loc=field%mesh%faces%shape%loc
    
  end function face_loc_tensor
  
  pure function face_ngi_vector(field, face_number) result (face_ngi)
    ! Return the number of nodes of element face_number.
    integer :: face_ngi
    type(vector_field),intent(in) :: field
    integer, intent(in) :: face_number
    
    face_ngi=field%mesh%faces%shape%ngi
    
  end function face_ngi_vector
  
  elemental function face_ele_vector(field, face_number) result (face_ele)
    ! Return the index of the element of which face is a part.
    integer :: face_ele
    type(vector_field), intent(in) :: field
    integer, intent(in) :: face_number

    face_ele=field%mesh%faces%face_element_list(face_number)

  end function face_ele_vector
  
  elemental function face_ele_mesh(mesh, face_number) result (face_ele)
    ! Return the index of the element of which face is a part.
    integer :: face_ele
    type(mesh_type), intent(in) :: mesh
    integer, intent(in) :: face_number

    face_ele=mesh%faces%face_element_list(face_number)

  end function face_ele_mesh
  
 function face_local_nodes_mesh(mesh, face_number) result (face_nodes)
    ! Return a pointer to a vector containing the local node numbers of
    ! facet face_number in mesh.
    integer, dimension(:), pointer :: face_nodes
    type(mesh_type),intent(in) :: mesh
    integer, intent(in) :: face_number
    
    ! This just reduces notational clutter.
    type(mesh_faces), pointer :: faces

    faces=>mesh%faces

    face_nodes=>faces%face_lno(faces%shape%loc*(face_number-1)+1:&
         &faces%shape%loc*face_number)
  
  end function face_local_nodes_mesh
  
  pure function surface_element_count_vector(field) result (element_count)
    ! Return the number of surface elements in a field.
    integer :: element_count
    type(vector_field),intent(in) :: field

    if (associated(field%mesh%faces)) then
      element_count=size(field%mesh%faces%boundary_ids)
    else
      element_count=0
    end if

  end function surface_element_count_vector
  
  pure function surface_element_count_mesh(mesh) result (element_count)
    ! Return the number of surface elements in a mesh.
    integer :: element_count
    type(mesh_type),intent(in) :: mesh

    if (associated(mesh%faces)) then
      element_count=size(mesh%faces%boundary_ids)
    else
      element_count=0
    end if

  end function surface_element_count_mesh  
  
  pure function face_count_scalar(field) result (face_count)
    ! Return the number of faces in a mesh.
    integer :: face_count
    type(scalar_field),intent(in) :: field

    if (associated(field%mesh%faces)) then
      face_count=size(field%mesh%faces%face_element_list)
    else
      face_count=0
    end if

  end function face_count_scalar
  
  pure function face_count_vector(field) result (face_count)
    ! Return the number of faces in a mesh.
    integer :: face_count
    type(vector_field),intent(in) :: field

    if (associated(field%mesh%faces)) then
      face_count=size(field%mesh%faces%face_element_list)
    else
      face_count=0
    end if

  end function face_count_vector
  
  pure function face_count_tensor(field) result (face_count)
    ! Return the number of faces in a mesh.
    integer :: face_count
    type(tensor_field),intent(in) :: field

    if (associated(field%mesh%faces)) then
      face_count=size(field%mesh%faces%face_element_list)
    else
      face_count=0
    end if

  end function face_count_tensor

  pure function face_count_mesh(mesh) result (face_count)
    ! Return the number of faces in a mesh.
    integer :: face_count
    type(mesh_type),intent(in) :: mesh

    if (associated(mesh%faces)) then
      face_count=size(mesh%faces%face_element_list)
    else
      face_count=0
    end if

  end function face_count_mesh  
  
  function face_local_nodes_vector(field, face_number) result (face_nodes)
    !!< Return a vector containing the local node numbers of
    !!< facet face_number in field.
    type(vector_field),intent(in) :: field
    integer, intent(in) :: face_number
    integer, dimension(face_loc(field, face_number)) :: face_nodes 

    face_nodes=face_local_nodes_mesh(field%mesh, face_number)

  end function face_local_nodes_vector
  
  function face_global_nodes_mesh(mesh, face_number) result (face_nodes)
    !!< Return a vector containing the global node numbers of
    !!< facet face_number in mesh.
    type(mesh_type),intent(in) :: mesh
    integer, intent(in) :: face_number
    integer, dimension(face_loc(mesh, face_number)) :: face_nodes 
    
    integer, dimension(:), pointer :: ele

    assert(has_faces(mesh))

    ele=>ele_nodes(mesh, face_ele(mesh,face_number))

    face_nodes=ele(face_local_nodes_mesh(mesh, face_number))
    
  end function face_global_nodes_mesh
    
  function face_global_nodes_scalar(field, face_number) result (face_nodes)
    !!< Return a vector containing the global node numbers of
    !!< facet face_number in field.
    type(scalar_field),intent(in) :: field
    integer, intent(in) :: face_number
    integer, dimension(face_loc(field, face_number)) :: face_nodes 

    face_nodes=face_global_nodes_mesh(field%mesh, face_number)

  end function face_global_nodes_scalar

  function face_global_nodes_vector(field, face_number) result (face_nodes)
    !!< Return a vector containing the global node numbers of
    !!< facet face_number in field.
    type(vector_field),intent(in) :: field
    integer, intent(in) :: face_number
    integer, dimension(face_loc(field, face_number)) :: face_nodes 

    face_nodes=face_global_nodes_mesh(field%mesh, face_number)

  end function face_global_nodes_vector

  function face_global_nodes_tensor(field, face_number) result (face_nodes)
    !!< Return a vector containing the global node numbers of
    !!< facet face_number in field.
    type(tensor_field),intent(in) :: field
    integer, intent(in) :: face_number
    integer, dimension(face_loc(field, face_number)) :: face_nodes 

    face_nodes=face_global_nodes_mesh(field%mesh, face_number)

  end function face_global_nodes_tensor
  
  function ele_shape_mesh(mesh, ele_number) result (ele_shape)
    ! Return a pointer to the shape of element ele_number.
    type(element_type), pointer :: ele_shape
    type(mesh_type),intent(in), target :: mesh
    integer, intent(in) :: ele_number
    
    ele_shape=>mesh%shape
    
  end function ele_shape_mesh
  
  function ele_shape_vector(field, ele_number) result (ele_shape)
    ! Return a pointer to the shape of element ele_number.
    type(element_type), pointer :: ele_shape
    type(vector_field),intent(in), target :: field
    integer, intent(in) :: ele_number
    
    ele_shape=>field%mesh%shape
    
  end function ele_shape_vector
  
  function ele_face_mesh(mesh, ele_number1, ele_number2) result (ele_face)
    ! Return the face in ele1 adjacent to ele2. 
    integer ele_face
    type(mesh_type), intent(in) :: mesh
    integer, intent(in) :: ele_number1, ele_number2

    ele_face=ival(mesh%faces%face_list, ele_number1, ele_number2)
    
  end function ele_face_mesh
  
  function ele_face_vector(field, ele_number1, ele_number2) result (ele_face)
    ! Return the face in ele1 adjacent to ele2. 
    integer ele_face
    type(vector_field), intent(in) :: field
    integer, intent(in) :: ele_number1, ele_number2
    
    ele_face=ele_face_mesh(field%mesh, ele_number1, ele_number2)

  end function ele_face_vector

  function ele_nodes_mesh(mesh, ele_number) result (ele_nodes)
    ! Return a pointer to a vector containing the global node numbers of
    ! element ele_number in mesh.
    integer, dimension(:), pointer :: ele_nodes
    type(mesh_type),intent(in) :: mesh
    integer, intent(in) :: ele_number

    ele_nodes=>mesh%ndglno(mesh%shape%loc*(ele_number-1)+1:&
         &mesh%shape%loc*ele_number)
  
  end function ele_nodes_mesh

  function ele_nodes_scalar(field, ele_number) result (ele_nodes)
    ! Return a pointer to a vector containing the global node numbers of
    ! element ele_number in field.
    integer, dimension(:), pointer :: ele_nodes
    type(scalar_field),intent(in) :: field
    integer, intent(in) :: ele_number

    ele_nodes=>ele_nodes_mesh(field%mesh, ele_number)

  end function ele_nodes_scalar

  function ele_nodes_vector(field, ele_number) result (ele_nodes)
    ! Return a pointer to a vector containing the global node numbers of
    ! element ele_number in field.
    integer, dimension(:), pointer :: ele_nodes
    type(vector_field),intent(in) :: field
    integer, intent(in) :: ele_number

    ele_nodes=>ele_nodes_mesh(field%mesh, ele_number)

  end function ele_nodes_vector

  function ele_nodes_tensor(field, ele_number) result (ele_nodes)
    ! Return a pointer to a tensor containing the global node numbers of
    ! element ele_number in field.
    integer, dimension(:), pointer :: ele_nodes
    type(tensor_field),intent(in) :: field
    integer, intent(in) :: ele_number

    ele_nodes=>ele_nodes_mesh(field%mesh, ele_number)

  end function ele_nodes_tensor

  pure function ele_loc_mesh(mesh, ele_number) result (ele_loc)
    ! Return the number of nodes of element ele_number.
    integer :: ele_loc
    type(mesh_type),intent(in) :: mesh
    integer, intent(in) :: ele_number
    
    ele_loc=mesh%shape%loc

  end function ele_loc_mesh

  pure function ele_loc_scalar(field, ele_number) result (ele_loc)
    ! Return the number of nodes of element ele_number.
    integer :: ele_loc
    type(scalar_field),intent(in) :: field
    integer, intent(in) :: ele_number
    
    ele_loc=field%mesh%shape%loc
    
  end function ele_loc_scalar

  pure function ele_loc_vector(field, ele_number) result (ele_loc)
    ! Return the number of nodes of element ele_number.
    integer :: ele_loc
    type(vector_field),intent(in) :: field
    integer, intent(in) :: ele_number
    
    ele_loc=field%mesh%shape%loc
    
  end function ele_loc_vector

  pure function ele_loc_tensor(field, ele_number) result (ele_loc)
    ! Return the number of nodes of element ele_number.
    integer :: ele_loc
    type(tensor_field),intent(in) :: field
    integer, intent(in) :: ele_number
    
    ele_loc=field%mesh%shape%loc
    
  end function ele_loc_tensor
  
  pure function face_vertices_mesh(mesh, face_number) result (face_vertices)
    ! Return the number of vertices of face face_number.
    integer :: face_vertices
    type(mesh_type),intent(in) :: mesh
    integer, intent(in) :: face_number
    
    face_vertices=mesh%faces%shape%quadrature%vertices

  end function face_vertices_mesh
  
  function face_vertices_shape(shape) result(vert)
    type(element_type), intent(in) :: shape
    integer :: vert

    select case(shape%numbering%type)
    case(ELEMENT_LAGRANGIAN, ELEMENT_BUBBLE, ELEMENT_TRACE)
      select case(shape%numbering%family)
      case (FAMILY_SIMPLEX)
        vert = shape%dim
      case (FAMILY_CUBE)
        vert = 2**(shape%dim-1)
      case default
        FLAbort("Unknown element family.")
      end select
    case default
      FLAbort("Unknown element type.")
    end select
    
  end function face_vertices_shape
  
  function face_shape_mesh(mesh, face_number) result (face_shape)
    ! Return a pointer to the shape of element ele_number.
    type(element_type), pointer :: face_shape
    type(mesh_type),intent(in) :: mesh
    integer, intent(in) :: face_number
    
    face_shape=>mesh%faces%shape
    
  end function face_shape_mesh

  function face_shape_scalar(field, face_number) result (face_shape)
    ! Return a pointer to the shape of element ele_number.
    type(element_type), pointer :: face_shape
    type(scalar_field),intent(in) :: field
    integer, intent(in) :: face_number
    
    face_shape=>field%mesh%faces%shape
    
  end function face_shape_scalar

  function face_shape_vector(field, face_number) result (face_shape)
    ! Return a pointer to the shape of element ele_number.
    type(element_type), pointer :: face_shape
    type(vector_field),intent(in) :: field
    integer, intent(in) :: face_number
    
    face_shape=>field%mesh%faces%shape
    
  end function face_shape_vector

  function face_shape_tensor(field, face_number) result (face_shape)
    ! Return a pointer to the shape of element ele_number.
    type(element_type), pointer :: face_shape
    type(tensor_field),intent(in) :: field
    integer, intent(in) :: face_number
    
    face_shape=>field%mesh%faces%shape
    
  end function face_shape_tensor
  
  function ele_val_scalar(field, ele_number) result (ele_val_out)
    ! Return the values of field at the nodes of ele_number.
    type(scalar_field),intent(in) :: field
    integer, intent(in) :: ele_number
    real, dimension(field%mesh%shape%loc) :: ele_val_out
    integer :: i

    select case(field%field_type)
    case(FIELD_TYPE_NORMAL)
      ele_val_out=field%val(ele_nodes(field,ele_number))
    case(FIELD_TYPE_CONSTANT)
      ele_val_out=field%val(1)
    case(FIELD_TYPE_PYTHON)
       call val_python
    end select
    
  contains
    
    subroutine val_python
      !!< This subroutine only exists to remove the following stack variables
      !!< from the main routine. 
      real, dimension(field%py_dim, field%mesh%shape%loc) :: pos
      real, dimension(field%py_dim, field%py_positions_shape%loc) :: tmp_pos

      if (.not. field%py_positions_same_mesh) then
        tmp_pos = ele_val(field%py_positions, ele_number)
        do i=1,field%py_dim
          pos(i, :) = matmul(field%py_locweight, tmp_pos(i, :))
        end do
      else
        do i=1,field%py_dim
          pos(i, :) = field%py_positions%val(i,ele_nodes(field%py_positions%mesh, ele_number))
        end do
      end if
      call set_from_python_function(ele_val_out, trim(field%py_func), pos, &
        time=current_time)

    end subroutine val_python

  end function ele_val_scalar

  function ele_val_vector(field, ele_number) result (ele_val)
    ! Return the values of field at the nodes of ele_number.
    type(vector_field),intent(in) :: field
    integer, intent(in) :: ele_number
    real, dimension(field%dim, field%mesh%shape%loc) :: ele_val

    integer :: i
    integer, dimension(:), pointer :: nodes
    
    select case(field%field_type)
    case(FIELD_TYPE_NORMAL)
      nodes => ele_nodes(field, ele_number)
      do i=1,field%dim
         ele_val(i, :) = field%val(i,nodes)
      end do
    case(FIELD_TYPE_CONSTANT)
      do i=1,field%dim
         ele_val(i,:)=field%val(i,1)
      end do
    end select


  end function ele_val_vector
  
  function ele_val_vector_dim(field, dim, ele_number) result (ele_val)
    ! Return the values of dimension dim of field at the nodes of ele_number.
    type(vector_field),intent(in) :: field
    integer, intent(in) :: ele_number
    real, dimension(field%mesh%shape%loc) :: ele_val
    integer, intent(in) :: dim

    select case(field%field_type)
    case(FIELD_TYPE_NORMAL)
      ele_val=field%val(dim,ele_nodes(field,ele_number))
    case(FIELD_TYPE_CONSTANT)
      ele_val=field%val(dim,1)
    end select

  end function ele_val_vector_dim

  function ele_val_tensor(field, ele_number) result (ele_val)
    ! Return the values of field at the nodes of ele_number.
    type(tensor_field),intent(in) :: field
    integer, intent(in) :: ele_number
    real, dimension(field%dim(1), field%dim(2), field%mesh%shape%loc) :: ele_val

    integer, dimension(:), pointer :: nodes
    integer :: i

    select case(field%field_type)
    case(FIELD_TYPE_NORMAL)
      nodes=>ele_nodes(field,ele_number)
      ele_val=field%val(:,:,nodes)
    case(FIELD_TYPE_CONSTANT)
      do i=1,size(ele_val, 3)
        ele_val(:, :, i)=field%val(:, :, 1)
      end do
    end select
    
  end function ele_val_tensor

  function ele_val_tensor_dim_dim(field, dim1, dim2, ele_number) result (ele_val)
    ! Return the values of field at the nodes of ele_number.
    type(tensor_field),intent(in) :: field
    integer, intent(in) :: dim1, dim2
    integer, intent(in) :: ele_number
    real, dimension(field%mesh%shape%loc) :: ele_val

    integer, dimension(:), pointer :: nodes

    select case(field%field_type)
    case(FIELD_TYPE_NORMAL)
      nodes=>ele_nodes(field,ele_number)
      ele_val=field%val(dim1,dim2,nodes)
    case(FIELD_TYPE_CONSTANT)
       ele_val=field%val(dim1, dim2, 1)
    end select
    
  end function ele_val_tensor_dim_dim
  
  function face_val_scalar(field, face_number) result (face_val)
    ! Return the values of field at the nodes of face_number.
    type(scalar_field),intent(in) :: field
    integer, intent(in) :: face_number
    real, dimension(face_loc(field, face_number)) :: face_val

    select case(field%field_type)
    case(FIELD_TYPE_NORMAL)
      face_val=field%val(face_global_nodes(field,face_number))
    case(FIELD_TYPE_CONSTANT)
      face_val=field%val(1)
    end select
    
  end function face_val_scalar

  function face_val_vector(field, face_number) result (face_val)
    ! Return the values of field at the nodes of face_number.
    type(vector_field),intent(in) :: field
    integer, intent(in) :: face_number
    real, dimension(field%dim, face_loc(field, face_number)) :: face_val

    integer :: i

    select case(field%field_type)
    case(FIELD_TYPE_NORMAL)
      face_val=field%val(:,face_global_nodes(field,face_number))
    case(FIELD_TYPE_CONSTANT)
      do i=1,field%dim
         face_val(i,:)=field%val(i,1)
      end do
    end select

  end function face_val_vector

  function face_val_vector_dim(field, dim, face_number) result (face_val)
    !!< Return the values of dimension dim of field at the nodes of face_number.
    type(vector_field),intent(in) :: field
    integer, intent(in) :: face_number
    real, dimension(face_loc(field, face_number)) :: face_val
    integer :: dim

    select case(field%field_type)
    case(FIELD_TYPE_NORMAL)
      face_val=field%val(dim,face_global_nodes(field,face_number))
    case(FIELD_TYPE_CONSTANT)
      face_val=field%val(dim,1)
    end select

  end function face_val_vector_dim

  function face_val_tensor(field, face_number) result (face_val)
    ! Return the values of field at the nodes of face_number.
    type(tensor_field),intent(in) :: field
    integer, intent(in) :: face_number
    real, dimension(field%dim(1), field%dim(2), face_loc(field, face_number)) ::&
         & face_val

    integer :: i

    select case(field%field_type)
    case(FIELD_TYPE_NORMAL)
       face_val=field%val(:,:,face_global_nodes(field,face_number))

    case(FIELD_TYPE_CONSTANT)
       forall(i=1:face_loc(field, face_number))
          face_val(:,:,i)=field%val(:,:,1)
       end forall
    end select

  end function face_val_tensor

  function face_val_tensor_dim_dim(field, dim1, dim2, face_number) result&
       & (face_val)
    !!< Return the values of dimension dim of field at the nodes of face_number.
    type(tensor_field),intent(in) :: field
    integer, intent(in) :: dim1, dim2
    integer, intent(in) :: face_number
    real, dimension(face_loc(field, face_number)) :: face_val

    select case(field%field_type)
    case(FIELD_TYPE_NORMAL)
      face_val=field%val(dim1,dim2,face_global_nodes(field,face_number))
    case(FIELD_TYPE_CONSTANT)
      face_val=field%val(dim1,dim2,1)
    end select

  end function face_val_tensor_dim_dim

  function node_val_scalar(field, node_number) result (val)
    ! Return the value of field at node node_number
    type(scalar_field),intent(in) :: field
    integer, intent(in) :: node_number
    real :: val

    select case(field%field_type)
    case(FIELD_TYPE_NORMAL)
      val=field%val(node_number)
    case(FIELD_TYPE_CONSTANT)
      val=field%val(1)
    case(FIELD_TYPE_PYTHON)
        FLAbort("node_val_scalar: FIELD_TYPE_PYTHON not implemented.")
!       call val_python
    end select

  end function node_val_scalar

  pure function node_val_vector(field, node_number) result (val)
    ! Return the value of field at node node_number
    type(vector_field),intent(in) :: field
    integer, intent(in) :: node_number
    real, dimension(field%dim) :: val

    integer :: i

    select case(field%field_type)
    case(FIELD_TYPE_NORMAL)
      do i=1,field%dim
         val(i)=field%val(i,node_number)
      end do
    case(FIELD_TYPE_CONSTANT)
      do i=1,field%dim
        val(i)=field%val(i,1)
      end do
    end select

  end function node_val_vector

  pure function node_val_tensor(field, node_number) result (val)
    ! Return the value of field at node node_number
    type(tensor_field),intent(in) :: field
    integer, intent(in) :: node_number
    real, dimension(field%dim(1), field%dim(2)) :: val

    select case(field%field_type)
    case(FIELD_TYPE_NORMAL)
      val=field%val(:,:,node_number)
    case(FIELD_TYPE_CONSTANT)
      val=field%val(:,:,1)
    end select

  end function node_val_tensor

  pure function node_val_tensor_dim_dim(field, dim1, dim2, node_number) result (val)
    ! Return the value of field at node node_number
    type(tensor_field),intent(in) :: field
    integer, intent(in) :: node_number
    integer, intent(in) :: dim1, dim2
    real :: val

    select case(field%field_type)
    case(FIELD_TYPE_NORMAL)
      val=field%val(dim1,dim2,node_number)
    case(FIELD_TYPE_CONSTANT)
      val=field%val(dim1,dim2,1)
    end select

  end function node_val_tensor_dim_dim

  pure function node_val_tensor_dim_dim_v(field, dim1, dim2, node_numbers) result (val)
    ! Return the value of field at nodes node_numbers
    type(tensor_field),intent(in) :: field
    integer, dimension(:), intent(in) :: node_numbers
    integer, intent(in) :: dim1, dim2
    real, dimension(size(node_numbers)) :: val

    select case(field%field_type)
    case(FIELD_TYPE_NORMAL)
      val=field%val(dim1,dim2,node_numbers)
    case(FIELD_TYPE_CONSTANT)
      val=field%val(dim1,dim2,1)
    end select

  end function node_val_tensor_dim_dim_v

  pure function node_val_scalar_v(field, node_numbers) result (val)
    ! Return the value of field at node node_numbers
    type(scalar_field),intent(in) :: field
    integer, dimension(:), intent(in) :: node_numbers
    real, dimension(size(node_numbers)) :: val

    select case(field%field_type)
    case(FIELD_TYPE_NORMAL)
      val=field%val(node_numbers)
    case(FIELD_TYPE_CONSTANT)
      val=field%val(1)
    end select

  end function node_val_scalar_v

  pure function node_val_vector_v(field, node_numbers) result (val)
    ! Return the value of field at node node_numbers
    type(vector_field),intent(in) :: field
    integer, dimension(:), intent(in) :: node_numbers
    real, dimension(field%dim, size(node_numbers)) :: val

    integer :: i

    select case(field%field_type)
    case(FIELD_TYPE_NORMAL)
      do i=1,field%dim
         val(i,:)=field%val(i,node_numbers)
      end do
    case(FIELD_TYPE_CONSTANT)
      do i=1,field%dim
         val(i,:)=field%val(i,1)
      end do
    end select

  end function node_val_vector_v

  pure function node_val_vector_dim(field, dim, node_number) &
       result (val) 
    ! Return the value of field at node node_numbers
    type(vector_field),intent(in) :: field
    integer, intent(in) :: node_number
    integer, intent(in) :: dim
    real :: val

    select case(field%field_type)
    case(FIELD_TYPE_NORMAL)
      val = field%val(dim,node_number)
    case(FIELD_TYPE_CONSTANT)
      val = field%val(dim,1)
    end select
    
  end function node_val_vector_dim

  pure function node_val_vector_dim_v(field, dim, node_numbers) &
       result (val) 
    ! Return the value of field at node node_numbers
    type(vector_field),intent(in) :: field
    integer, dimension(:), intent(in) :: node_numbers
    integer, intent(in) :: dim
    real, dimension(size(node_numbers)) :: val

    select case(field%field_type)
    case(FIELD_TYPE_NORMAL)
      val(:)=field%val(dim,node_numbers)
    case(FIELD_TYPE_CONSTANT)
      val=field%val(dim,1)
    end select

  end function node_val_vector_dim_v

  pure function node_val_tensor_v(field, node_numbers) result (val)
    ! Return the value of field at node node_numbers
    type(tensor_field),intent(in) :: field
    integer, dimension(:), intent(in) :: node_numbers
    real, dimension(field%dim(1), field%dim(2), size(node_numbers)) :: val
    integer :: i


    select case(field%field_type)
    case(FIELD_TYPE_NORMAL)
      val=field%val(:,:,node_numbers)
    case(FIELD_TYPE_CONSTANT)
      do i=1,size(node_numbers)
        val(:, :, i)=field%val(:, :, 1)
      end do
    end select

  end function node_val_tensor_v
  
  ! ------------------------------------------------------------------------
  ! Point wise python evalutiation
  ! ------------------------------------------------------------------------
  ! these could go in a separate module
  
  subroutine set_values_from_python_scalar(values, func, x, y, z, time)
    !!< Given a list of positions and a time, evaluate the python function
    !!< specified in the string func at those points. 
    real, dimension(:), intent(inout) :: values
    !! Func may contain any python at all but the following function must
    !! be defiled:
    !!  def val(X, t)
    !! where X is a tuple containing the position of a point and t is the
    !! time. The result must be a float. 
    character(len=*), intent(in) :: func
    real, dimension(size(values)), target :: x
    real, dimension(size(values)), optional, target :: y
    real, dimension(size(values)), optional, target :: z
    real :: time

    real, dimension(:), pointer :: lx, ly, lz
    real, dimension(0), target :: zero
    integer :: stat, dim
    
    lx=>x
    ly=>zero
    lz=>zero
    dim=1
    if (present(y)) then
       ly=>y
       dim=2
       if (present(z)) then
          lz=>z
          dim=3
       end if
    end if
    
    call set_scalar_field_from_python(func, len(func), dim,&
            & size(values), lx, ly, lz, time, values, stat)

    if (stat/=0) then
      ewrite(-1, *) "Python error, Python string was:"
      ewrite(-1, *) trim(func)
      FLAbort("Dying")
    end if

  end subroutine set_values_from_python_scalar
  
  subroutine set_values_from_python_scalar_pos(values, func, pos, time)
    !!< Given a list of positions and a time, evaluate the python function
    !!< specified in the string func at those points. 
    real, dimension(:), intent(inout) :: values
    !! Func may contain any python at all but the following function must
    !! be defiled:
    !!  def val(X, t)
    !! where X is a tuple containing the position of a point and t is the
    !! time. The result must be a float. 
    character(len=*), intent(in) :: func
    real, dimension(:, :), intent(in), target :: pos
    real :: time

    real, dimension(:), pointer :: lx, ly, lz
    real, dimension(0), target :: zero
    integer :: stat, dim
    
    dim=size(pos, 1)
    select case(dim)
    case(1)
      lx=>pos(1, :)
      ly=>zero
      lz=>zero
    case(2)
      lx=>pos(1, :)
      ly=>pos(2, :)
      lz=>zero
    case(3)
      lx=>pos(1, :)
      ly=>pos(2, :)
      lz=>pos(3, :)
    end select
    
    call set_scalar_field_from_python(func, len(func), dim,&
            & size(values), lx, ly, lz, time, values, stat)

    if (stat/=0) then
      ewrite(-1, *) "Python error, Python string was:"
      ewrite(-1, *) trim(func)
      FLAbort("Dying")
    end if
  end subroutine set_values_from_python_scalar_pos

  subroutine set_values_from_python_vector(values, func, x, y, z, time)
    !!< Given a list of positions and a time, evaluate the python function
    !!< specified in the string func at those points. 
    real, dimension(:,:), target, intent(inout) :: values
    !! Func may contain any python at all but the following function must
    !! be defiled:
    !!  def val(X, t)
    !! where X is a tuple containing the position of a point and t is the
    !! time. The result must be a float. 
    character(len=*), intent(in) :: func
    real, dimension(size(values,2)), target :: x
    real, dimension(size(values,2)), optional, target :: y
    real, dimension(size(values,2)), optional, target :: z
    real :: time

    real, dimension(:), pointer :: lx, ly, lz
    real, dimension(:), pointer :: lvx,lvy,lvz
    real, dimension(0), target :: zero
    integer :: stat, dim
    
    lx=>x
    ly=>zero
    lz=>zero
    dim=1
    if (present(y)) then
       ly=>y
       dim=2
       if (present(z)) then
          lz=>z
          dim=3
       end if
    end if

    lvx=>values(1,:)
    lvy=>zero
    lvz=>zero
    if(size(values,1)>1) then
       lvy=>values(2,:)
       if(size(values,1)>2) then
          lvz => values(3,:)
       end if
    end if
    call set_vector_field_from_python(func, len_trim(func), dim,&
            & size(values,2), lx, ly, lz, time,size(values,1), &
            lvx,lvy,lvz, stat)

    if (stat/=0) then
      ewrite(-1, *) "Python error, Python string was:"
      ewrite(-1, *) trim(func)
      FLAbort("Dying")
    end if

  end subroutine set_values_from_python_vector

  subroutine set_values_from_python_vector_pos(values, func, pos, time)
    !!< Given a list of positions and a time, evaluate the python function
    !!< specified in the string func at those points. 
    real, dimension(:,:), intent(inout) :: values
    !! Func may contain any python at all but the following function must
    !! be defiled:
    !!  def val(X, t)
    !! where X is a tuple containing the position of a point and t is the
    !! time. The result must be a float. 
    character(len=*), intent(in) :: func
    real, dimension(:, :), intent(in), target :: pos
    real, intent(in) :: time

    integer :: dim
    
    dim=size(pos, 1)
    select case(dim)
    case(1)
      call set_values_from_python_vector(values, func, pos(1,:), time=time)
    case(2)
      call set_values_from_python_vector(values, func, pos(1,:), pos(2,:), time=time)
    case(3)
      call set_values_from_python_vector(values, func, pos(1,:), pos(2,:), pos(3,:), time=time)
    end select    
    
  end subroutine set_values_from_python_vector_pos

  subroutine set_values_from_python_vector_field(values, func, vfield, time)
    !!< Given a list of positions and a time, evaluate the python function
    !!< specified in the string func at those points. 
    real, dimension(:,:), intent(inout) :: values
    !! Func may contain any python at all but the following function must
    !! be defiled:
    !!  def val(X, t)
    !! where X is a tuple containing the position of a point and t is the
    !! time. The result must be a float. 
    character(len=*), intent(in) :: func
    type(vector_field), intent(in) :: vfield
    real, intent(in) :: time

    integer :: dim
    
    dim=vfield%dim
    select case(dim)
    case(1)
      call set_values_from_python_vector(values, func, vfield%val(1,:), time=time)
    case(2)
      call set_values_from_python_vector(values, func, vfield%val(1,:), vfield%val(2,:), time=time)
    case(3)
      call set_values_from_python_vector(values, func, vfield%val(1,:), vfield%val(2,:), vfield%val(3,:), time=time)
    end select    
    
  end subroutine set_values_from_python_vector_field

  function local_face_number_mesh(mesh, global_face_number, stat) result (local_face_number)
    ! Return the local face number of the given global face number in element ele_number
    integer :: local_face_number
    type(mesh_type), intent(in) :: mesh
    integer, intent(in) :: global_face_number
    integer, intent(inout), optional :: stat
    
    integer :: ele_number, i
    integer, dimension(:), pointer :: element_faces
    
    local_face_number = 0
    
    ele_number = face_ele(mesh, global_face_number)
    
    element_faces => ele_faces(mesh, ele_number)
    do i = 1, size(element_faces)
      if(element_faces(i) == global_face_number) then
        local_face_number = i
        exit
      end if
    end do
    
    if(local_face_number==0) then
      if(present(stat)) then
        stat = 1
      else
        FLAbort("Failed to find local face number.")
      end if
    else
      if(present(stat)) stat = 0
    end if
    
  end function local_face_number_mesh
  
  pure function element_count_scalar(field) result (element_count)
    ! Return the number of elements in a field.
    integer :: element_count
    type(scalar_field),intent(in) :: field

    element_count=field%mesh%elements

  end function element_count_scalar

  pure function element_count_vector(field) result (element_count)
    ! Return the number of elements in a field.
    integer :: element_count
    type(vector_field),intent(in) :: field

    element_count=field%mesh%elements

  end function element_count_vector

  pure function element_count_tensor(field) result (element_count)
    ! Return the number of elements in a field.
    integer :: element_count
    type(tensor_field),intent(in) :: field

    element_count=field%mesh%elements

  end function element_count_tensor
  
  pure function element_count_mesh(mesh) result (element_count)
    ! Return the number of nodes in a mesh.
    integer :: element_count
    type(mesh_type),intent(in) :: mesh

    element_count=mesh%elements

  end function element_count_mesh  
  
  pure function node_count_mesh(mesh) result (node_count)
    ! Return the number of nodes in a mesh.
    integer :: node_count
    type(mesh_type), intent(in) :: mesh

    node_count=mesh%nodes

  end function node_count_mesh

  pure function node_count_scalar(field) result (node_count)
    ! Return the number of nodes in a field.
    integer :: node_count
    type(scalar_field),intent(in) :: field

    node_count=node_count_mesh(field%mesh)

  end function node_count_scalar

  pure function node_count_vector(field) result (node_count)
    ! Return the number of nodes in a field.
    integer :: node_count
    type(vector_field),intent(in) :: field

    node_count=node_count_mesh(field%mesh)

  end function node_count_vector

  pure function node_count_tensor(field) result (node_count)
    ! Return the number of nodes in a field.
    integer :: node_count
    type(tensor_field),intent(in) :: field

    node_count=node_count_mesh(field%mesh)

  end function node_count_tensor
  
  function has_faces_mesh(mesh) result (has_faces)
    ! Check whether the faces component of mesh has been calculated.
    logical :: has_faces
    type(mesh_type), intent(in) :: mesh

    has_faces=associated(mesh%faces)

  end function has_faces_mesh
  
  function extract_scalar_field_from_vector_field(vfield, dim, stat) result(sfield)
  !!< This function gives you a way to treat a vector field
  !!< as a union of scalar fields.

    type(vector_field), intent(in), target :: vfield
    integer, intent(in) :: dim
    type(scalar_field) :: sfield
    integer, intent(out), optional :: stat

    if (present(stat)) then
      stat = 0
      if (dim > vfield%dim) then
        stat = 1
        return
      end if
    end if
    assert(dim .le. vfield%dim)

    ! Note that the reference count is not incremented as this is a
    ! borrowed field reference.
    sfield%mesh = vfield%mesh
    sfield%val  => vfield%val(dim,:)
    sfield%val_stride = vfield%dim
    sfield%option_path = vfield%option_path
    sfield%field_type = vfield%field_type
    write(sfield%name, '(a, i0)') trim(vfield%name) // "%", dim
    
    ! FIXME: make these the same as the vector field
    sfield%py_dim = mesh_dim(vfield%mesh)
    sfield%py_positions_shape => vfield%mesh%shape
    
    sfield%refcount => vfield%refcount

  end function extract_scalar_field_from_vector_field
  
  ! ------------------------------------------------------------------------
  ! Geometric element volume routines. These really ought to go somewhere
  ! else but tend to cause dependency problems when they do.
  ! ------------------------------------------------------------------------

  function tetvol_new(positions, ele) result(t)
    real :: t
    type(vector_field), intent(in) :: positions
    integer, intent(in) :: ele
    real, dimension(positions%dim, ele_loc(positions, ele)) :: pos

    pos = ele_val(positions, ele)
    t = tetvol_old(pos(1, :), pos(2, :), pos(3, :))
    
  end function tetvol_new
  
  real function tetvol_old( x, y, z )

    real x(4), y(4), z(4)
    real vol, x12, x13, x14, y12, y13, y14, z12, z13, z14
    !
    x12 = x(2) - x(1)
    x13 = x(3) - x(1)
    x14 = x(4) - x(1)
    y12 = y(2) - y(1)
    y13 = y(3) - y(1)
    y14 = y(4) - y(1)
    z12 = z(2) - z(1)
    z13 = z(3) - z(1)
    z14 = z(4) - z(1)
    !
    vol = x12*( y13*z14 - y14*z13 )  &
         + x13*( y14*z12 - y12*z14 ) &
         + x14*( y12*z13 - y13*z12 )
    !
    tetvol_old = vol/6
    !
    return
  end function tetvol_old
  
  function triarea(positions, ele) result(t)
    type(vector_field), intent(in) :: positions
    integer, intent(in) :: ele
    real :: t
    real, dimension(positions%dim, positions%mesh%shape%loc) :: pos
    real :: xA, xB, yA, yB, xC, yC

    pos = ele_val(positions, ele)
    if (positions%dim == 2) then
      xA = pos(1, 1); xB = pos(1, 2); xC = pos(1, 3)
      yA = pos(2, 1); yB = pos(2, 2); yC = pos(2, 3)
      t = abs((xB*yA-xA*yB)+(xC*yB-xB*yC)+(xA*yC-xC*yA))/2
    elseif (positions%dim == 3) then
      ! http://mathworld.wolfram.com/TriangleArea.html, (19)
      t = 0.5 * norm2(cross_product(pos(:, 2) - pos(:, 1), pos(:, 1) - pos(:, 3)))
    else
      FLAbort("Only 2 or 3 dimensions supported, sorry")
    end if
  end function triarea
  
  function tetvol_1d(positions, ele) result(t)
    type(vector_field), intent(in) :: positions
    integer, intent(in) :: ele
    
    real :: t
    
    integer, dimension(:), pointer :: element_nodes => null()
    
    assert(positions%dim == 1)
    
    element_nodes => ele_nodes(positions, ele)
    
    assert(size(element_nodes) == 2)
    t = abs(node_val(positions, 1, element_nodes(2)) - node_val(positions, 1, element_nodes(1)))
  
  end function tetvol_1d
  
  function simplex_volume(positions, ele) result(t)
    type(vector_field), intent(in) :: positions
    integer, intent(in) :: ele
    real :: t

    select case(mesh_dim(positions))
      case(3)
        t = tetvol_new(positions, ele)
      case(2)
        t = triarea(positions, ele)
      case(1)
        t = tetvol_1d(positions, ele)
    case default
      FLAbort("Invalid dimension")
    end select
   
  end function simplex_volume

end module libsupermesh_fields_base
