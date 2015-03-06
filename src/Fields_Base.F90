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
!  use shape_functions, only: element_type	! IAKOVOS commented out
!  use tensors					! IAKOVOS commented out
  use libsupermesh_fields_data_types
  use libsupermesh_reference_counting
  use libsupermesh_global_parameters, only: FIELD_NAME_LEN, current_debug_level, current_time
  use libsupermesh_elements
  use libsupermesh_element_numbering
!  use embed_python				! IAKOVOS commented out
  use libsupermesh_sparse_tools
!  use vector_tools, only: solve, invert, norm2, cross_product	! IAKOVOS commented out
  implicit none
  
  interface ele_nodes
     module procedure ele_nodes_scalar, ele_nodes_vector, ele_nodes_tensor,&
          & ele_nodes_mesh
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
  
  interface ele_val
     module procedure ele_val_scalar, ele_val_vector, ele_val_vector_dim,&
          & ele_val_tensor, ele_val_tensor_dim_dim
  end interface
  
  interface ele_count
     module procedure element_count_scalar, element_count_vector,&
          & element_count_tensor, element_count_mesh
  end interface
  
  interface node_count
     module procedure node_count_scalar, node_count_vector,&
          & node_count_tensor, node_count_mesh
  end interface
  
  interface node_val
!     module procedure node_val_scalar, node_val_vector, node_val_tensor, &
!          & node_val_scalar_v, node_val_vector_v, node_val_vector_dim_v,&
!          & node_val_tensor_v, node_val_vector_dim, node_val_tensor_dim_dim, &
!          & node_val_tensor_dim_dim_v
     module procedure node_val_vector
  end interface
  
contains

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

end module libsupermesh_fields_base
