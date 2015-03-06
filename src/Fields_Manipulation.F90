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
module libsupermesh_fields_manipulation
use libsupermesh_elements
!use element_set		! IAKOVOS commented out
!use embed_python		! IAKOVOS commented out
use libsupermesh_data_structures
use libsupermesh_fields_data_types
use libsupermesh_fields_base
use libsupermesh_fields_allocates
use libsupermesh_halo_data_types
use libsupermesh_halos_allocates
use libsupermesh_halos_base
!use halos_debug		! IAKOVOS commented out
!use halos_numbering		! IAKOVOS commented out
!use halos_ownership		! IAKOVOS commented out
!use halos_repair		! IAKOVOS commented out
use libsupermesh_quicksort
use libsupermesh_parallel_tools
!use vector_tools		! IAKOVOS commented out
!use memory_diagnostics		! IAKOVOS commented out
implicit none

  private
  
  public :: piecewise_constant_mesh, set
  
  interface set
    module procedure set_vector_field_node, set_vector_field, &
                   & set_vector_field_node_dim, set_vector_field_dim, &
                   & set_vector_field_nodes, &
                   & set_vector_field_nodes_dim, &
                   & set_vector_field_field, &
                   & set_vector_field_field_dim, &
                   & set_vector_field_theta, &
                   & set_vector_field_vfield_dim, &
                   & set_scalar_field_node
  end interface
    
  contains

  function piecewise_constant_mesh(in_mesh, name) result(new_mesh)
    !!< From a given mesh, return a scalar field
    !!< allocated on the mesh that's topologically the same
    !!< but has piecewise constant basis functions.
    !!< This is for the definition of elementwise quantities.
    type(mesh_type), intent(in) :: in_mesh
    type(mesh_type) :: new_mesh
    type(element_type) :: shape, old_shape
    character(len=*), intent(in) :: name

    old_shape = in_mesh%shape

!    shape = make_element_shape(vertices=old_shape%numbering%vertices, dim=old_shape%dim, degree=0, quad=old_shape%quadrature)		! ToDo 1
!    new_mesh = make_mesh(model=in_mesh, shape=shape, continuity=-1)									! ToDo 1
    new_mesh%name=name
!    call deallocate(shape)														! ToDo 1
    
  end function piecewise_constant_mesh
  
  subroutine set_vector_field_node(field, node, val)
    !!< Set the vector field at the specified node
    !!< Does not work for constant fields
    type(vector_field), intent(inout) :: field
    integer, intent(in) :: node
    real, intent(in), dimension(:) :: val
    integer :: i

    assert(field%field_type==FIELD_TYPE_NORMAL)
    
    do i=1,field%dim
      field%val(i,node) = val(i)
    end do
    
  end subroutine set_vector_field_node
  
  subroutine set_vector_field_node_dim(field, dim, node, val)
    !!< Set the vector field at the specified node
    !!< Does not work for constant fields
    type(vector_field), intent(inout) :: field
    integer, intent(in) :: node
    real, intent(in) :: val
    integer, intent(in):: dim

    assert(field%field_type==FIELD_TYPE_NORMAL)
    assert(dim>=1 .and. dim<=field%dim)
    
    field%val(dim,node) = val
    
  end subroutine set_vector_field_node_dim
  
  subroutine set_vector_field_nodes(field, node_numbers, val)
    !!< Set the vector field at the specified nodes
    !!< Does not work for constant fields
    type(vector_field), intent(inout) :: field
    integer, dimension(:), intent(in) :: node_numbers
    !! values to set ( dimension x #nodes)
    real, intent(in), dimension(:,:) :: val
    integer :: i

    assert(field%field_type==FIELD_TYPE_NORMAL)
    
    do i=1,field%dim
      field%val(i,node_numbers) = val(i, :)
    end do
    
  end subroutine set_vector_field_nodes
  
  subroutine set_vector_field_nodes_dim(field, dim, node_numbers, val)
    !!< Set the vector field at the specified nodes
    !!< Does not work for constant fields
    type(vector_field), intent(inout) :: field
    integer, dimension(:), intent(in) :: node_numbers
    !! values to set
    real, intent(in), dimension(:) :: val
    integer, intent(in) :: dim

    assert(field%field_type==FIELD_TYPE_NORMAL)
    assert(dim>=1 .and. dim<=field%dim)
    
    field%val(dim,node_numbers) = val
    
  end subroutine set_vector_field_nodes_dim
  
  subroutine set_vector_field(field, val)
    !!< Set the vector field with a constant value
    !!< Works for constant and space varying fields.
    type(vector_field), intent(inout) :: field
    real, intent(in), dimension(:) :: val
    integer :: i

    assert(field%field_type/=FIELD_TYPE_PYTHON)
    
    do i=1,field%dim
      field%val(i,:) = val(i)
    end do

  end subroutine set_vector_field
  
  subroutine set_vector_field_dim(field, dim, val)
    !!< Set the vector field with a constant value
    !!< Works for constant and space varying fields.
    type(vector_field), intent(inout) :: field
    real, intent(in):: val
    integer, intent(in):: dim

    assert(field%field_type/=FIELD_TYPE_PYTHON)
    assert(dim>=1 .and. dim<=field%dim)
    
    field%val(dim,:) = val

  end subroutine set_vector_field_dim

  subroutine set_vector_field_arr(field, val)
    !!< Set the vector field with an array for all nodes at once
    type(vector_field), intent(inout) :: field
    real, intent(in), dimension(:, :) :: val
    integer :: i

    assert(field%field_type == FIELD_TYPE_NORMAL)
    
    do i=1,field%dim
      field%val(i,:) = val(i, :)
    end do

  end subroutine set_vector_field_arr

  subroutine set_vector_field_arr_dim(field, dim, val)
    !!< Set the vector field with an array for all nodes at once
    type(vector_field), intent(inout) :: field
    real, intent(in), dimension(:) :: val
    integer, intent(in):: dim
    
    assert(field%field_type == FIELD_TYPE_NORMAL)
    assert(dim>=1 .and. dim<=field%dim)
    
    field%val(dim,:) = val

  end subroutine set_vector_field_arr_dim

  subroutine set_vector_field_field(out_field, in_field )
    !!< Set in_field to out_field. This will only work if the fields have
    !!< the same mesh.
    type(vector_field), intent(inout) :: out_field
    type(vector_field), intent(in) :: in_field
    
    integer :: dim

    assert(mesh_compatible(out_field%mesh, in_field%mesh))
    assert(out_field%field_type/=FIELD_TYPE_PYTHON)
    assert(out_field%field_type==FIELD_TYPE_NORMAL .or. in_field%field_type==FIELD_TYPE_CONSTANT)
    assert(in_field%dim==out_field%dim)

    select case (in_field%field_type)
    case (FIELD_TYPE_NORMAL)
       do dim=1,in_field%dim
          out_field%val(dim,:)=in_field%val(dim,:)
       end do
    case (FIELD_TYPE_CONSTANT)
       do dim=1,in_field%dim
          out_field%val(dim,:)=in_field%val(dim,1)
       end do
    case default
       ! someone could implement in_field type python
!       FLAbort("Illegal in_field field type in set()")			! ToDo
    end select

  end subroutine set_vector_field_field

  subroutine set_vector_field_theta(out_field, in_field_new, in_field_old, theta)
    !!< Set theta*in_field_new + (1.-theta)*in_field_old to out_field. This will only work if the fields have
    !!< the same mesh.
    type(vector_field), intent(inout) :: out_field
    type(vector_field), intent(in) :: in_field_new, in_field_old
    real, intent(in) :: theta
    
    integer :: dim
    
    assert(mesh_compatible(out_field%mesh, in_field_new%mesh))
    assert(mesh_compatible(out_field%mesh, in_field_old%mesh))
    assert(out_field%field_type/=FIELD_TYPE_PYTHON)
#ifndef NDEBUG
    if(.not.(out_field%field_type==FIELD_TYPE_NORMAL .or. &
       (in_field_new%field_type==FIELD_TYPE_CONSTANT .and. &
        in_field_old%field_type==FIELD_TYPE_CONSTANT))) then
       ewrite(-1,*) "Incompatible field types in set()"
!        FLAbort("Evilness");						! ToDo
    end if
#endif
    assert(in_field_new%dim==out_field%dim)
    assert(in_field_old%dim==out_field%dim)
    
    select case (in_field_new%field_type)
    case (FIELD_TYPE_NORMAL)
      do dim = 1, out_field%dim
        out_field%val(dim,:)=theta*in_field_new%val(dim,:) + (1.-theta)*in_field_old%val(dim,:)
      end do
    case (FIELD_TYPE_CONSTANT)
      do dim = 1, out_field%dim
        out_field%val(dim,:)=theta*in_field_new%val(dim,1) + (1.-theta)*in_field_old%val(dim,1)
      end do
    case default
       ! someone could implement in_field type python
!       FLAbort("Illegal in_field field type in set()")			! ToDo
    end select

  end subroutine set_vector_field_theta

  subroutine set_vector_field_field_dim(out_field, dim, in_field)
    !!< Set in_field to out_field. This will only work if the fields have
    !!< the same mesh.
    type(vector_field), intent(inout) :: out_field
    type(scalar_field), intent(in) :: in_field
    integer, intent(in):: dim

    assert(mesh_compatible(out_field%mesh, in_field%mesh))
    assert(out_field%field_type/=FIELD_TYPE_PYTHON)
    assert(out_field%field_type==FIELD_TYPE_NORMAL.or.in_field%field_type==FIELD_TYPE_CONSTANT)
    assert(dim>=1 .and. dim<=out_field%dim)

    select case (in_field%field_type)
    case (FIELD_TYPE_NORMAL)
       out_field%val(dim,:)=in_field%val
    case (FIELD_TYPE_CONSTANT)
       out_field%val(dim,:)=in_field%val(1)
    case default
       ! someone could implement in_field type python
!       FLAbort("Illegal in_field field type in set()")			! ToDo
    end select

  end subroutine set_vector_field_field_dim

  subroutine set_vector_field_vfield_dim(out_field, dim, in_field)
    !!< Set in_field to out_field. This will only work if the fields have
    !!< the same mesh.
    type(vector_field), intent(inout) :: out_field
    type(vector_field), intent(in) :: in_field
    integer, intent(in):: dim

    assert(mesh_compatible(out_field%mesh, in_field%mesh))
    assert(out_field%field_type/=FIELD_TYPE_PYTHON)
    assert(out_field%field_type==FIELD_TYPE_NORMAL.or.in_field%field_type==FIELD_TYPE_CONSTANT)
    assert(dim>=1 .and. dim<=out_field%dim .and. dim<=in_field%dim)

    select case (in_field%field_type)
    case (FIELD_TYPE_NORMAL)
       out_field%val(dim,:)=in_field%val(dim,:)
    case (FIELD_TYPE_CONSTANT)
       out_field%val(dim,:)=in_field%val(dim,1)
    case default
       ! someone could implement in_field type python
!       FLAbort("Illegal in_field field type in set()")			! ToDo
    end select

  end subroutine set_vector_field_vfield_dim
  
  subroutine set_scalar_field_node(field, node_number, val)
    !!< Set the scalar field at the specified node
    !!< Does not work for constant fields
    type(scalar_field), intent(inout) :: field
    integer, intent(in) :: node_number
    real, intent(in) :: val

    assert(field%field_type==FIELD_TYPE_NORMAL)
    
    field%val(node_number) = val
    
  end subroutine set_scalar_field_node
  
end module libsupermesh_fields_manipulation
