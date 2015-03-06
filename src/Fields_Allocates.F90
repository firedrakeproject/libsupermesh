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
module libsupermesh_fields_allocates
use libsupermesh_elements
use libsupermesh_fields_data_types
use libsupermesh_fields_base
!use shape_functions, only: make_element_shape		! IAKOVOS commented out
use libsupermesh_global_parameters, only: PYTHON_FUNC_LEN, empty_path, empty_name, &
     topology_mesh_name, NUM_COLOURINGS
use libsupermesh_halo_data_types
use libsupermesh_halos_allocates
!use halos_repair		! IAKOVOS commented out
!use pickers_deallocates	! IAKOVOS commented out
!use adjacency_lists		! IAKOVOS commented out
!use global_numbering, only: make_global_numbering, make_global_numbering_dg,&
!     &make_global_numbering_trace	! IAKOVOS commented out
!use memory_diagnostics		! IAKOVOS commented out
!use ieee_arithmetic		! IAKOVOS commented out
use libsupermesh_data_structures
use libsupermesh_parallel_tools

implicit none

  private

! IAKOVOS commented out
!  public :: allocate, deallocate, incref, decref, has_references, add_faces, &
!     & deallocate_faces, zero
  public :: allocate, deallocate, incref, decref, has_references, &
    & zero
! IAKOVOS commented out
!  public :: make_element_shape, make_mesh, make_mesh_periodic, make_submesh, &
!    & create_surface_mesh, make_fake_mesh_linearnonconforming
! IAKOVOS commented out
!  public :: extract_scalar_field, wrap_mesh, wrap_scalar_field, &
!    & wrap_tensor_field
  public :: wrap_mesh
! IAKOVOS commented out
!  public :: add_lists, extract_lists, add_nnlist, extract_nnlist, add_nelist, &
!    & extract_nelist, add_eelist, extract_eelist, remove_lists, remove_nnlist, &
!    & remove_nelist, remove_eelist, extract_elements, remove_boundary_conditions
  public :: extract_eelist, remove_lists

  interface allocate
! IAKOVOS commented out
!     module procedure allocate_scalar_field, allocate_vector_field,&
!          & allocate_tensor_field, allocate_mesh, &
!          & allocate_scalar_boundary_condition, &
!          & allocate_vector_boundary_condition
     module procedure allocate_mesh
  end interface

! IAKOVOS commented out
!  interface deallocate
!     module procedure deallocate_mesh, deallocate_scalar_field,&
!          & deallocate_vector_field, deallocate_tensor_field, &
!          & deallocate_scalar_boundary_condition, &
!          & deallocate_vector_boundary_condition
!  end interface

  interface zero
     module procedure zero_scalar, zero_vector, zero_tensor, &
          zero_vector_dim, zero_tensor_dim_dim, &
          zero_scalar_field_nodes, zero_vector_field_nodes, zero_tensor_field_nodes
  end interface

! IAKOVOS commented out
!  interface deallocate_faces
!     module procedure deallocate_mesh_faces
!  end interface

! IAKOVOS commented out
!  interface add_lists
!    module procedure add_lists_mesh, add_lists_scalar, add_lists_vector, &
!      & add_lists_tensor
!  end interface add_lists

! IAKOVOS commented out
!  interface extract_lists
!    module procedure extract_lists_mesh, extract_lists_scalar, &
!      & extract_lists_vector, extract_lists_tensor
!  end interface extract_lists
 
! IAKOVOS commented out
!  interface add_nnlist
!    module procedure add_nnlist_mesh, add_nnlist_scalar, add_nnlist_vector, &
!      & add_nnlist_tensor
!  end interface add_nnlist

! IAKOVOS commented out  
!  interface extract_nnlist
!    module procedure extract_nnlist_mesh, extract_nnlist_scalar, &
!      & extract_nnlist_vector, extract_nnlist_tensor
!  end interface extract_nnlist

! IAKOVOS commented out
!  interface add_nelist
!    module procedure add_nelist_mesh, add_nelist_scalar, add_nelist_vector, &
!      & add_nelist_tensor
!  end interface add_nelist
 
! IAKOVOS commented out
!  interface extract_nelist
!    module procedure extract_nelist_mesh, extract_nelist_scalar, &
!      & extract_nelist_vector, extract_nelist_tensor
!  end interface extract_nelist
 
! IAKOVOS commented out
!  interface add_eelist
!    module procedure add_eelist_mesh, add_eelist_scalar, add_eelist_vector, &
!      & add_eelist_tensor
!  end interface add_eelist

! IAKOVOS commented out
  interface extract_eelist
!    module procedure extract_eelist_mesh, extract_eelist_scalar, &
!      & extract_eelist_vector, extract_eelist_tensor
    module procedure extract_eelist_vector, extract_eelist_mesh
  end interface extract_eelist
  
  interface remove_lists
    module procedure remove_lists_mesh
  end interface remove_lists
  
  ! IAKOVOS commented out
!  interface remove_nnlist
!    module procedure remove_nnlist_mesh
!  end interface remove_nnlist
  
! IAKOVOS commented out
!  interface remove_nelist
!    module procedure remove_nelist_mesh
!  end interface remove_nelist
  
  ! IAKOVOS commented out
!  interface remove_eelist
!    module procedure remove_eelist_mesh
!  end interface remove_eelist
 
! IAKOVOS commented out
!  interface remove_boundary_conditions
!    module procedure remove_boundary_conditions_scalar, &
!      remove_boundary_conditions_vector
!  end interface remove_boundary_conditions
  
#include "Reference_count_interface_mesh_type.F90"
#include "Reference_count_interface_scalar_field.F90"
#include "Reference_count_interface_vector_field.F90"
#include "Reference_count_interface_tensor_field.F90"

contains

  subroutine allocate_mesh(mesh, nodes, elements, shape, name)
    type(mesh_type), intent(out) :: mesh
    integer, intent(in) :: nodes, elements
    type(element_type), target, intent(in) :: shape
    character(len=*), intent(in), optional :: name
    integer :: i
#ifdef _OPENMP
    integer :: j
#endif
    
    mesh%nodes=nodes

    mesh%elements=elements

    mesh%shape=shape
    call incref(shape)
    
    if (present(name)) then
       mesh%name=name
    else
       mesh%name=empty_name
    end if
    
    ! should happen in derived type initialisation already,
    ! but just to make sure in case an mesh variable is supplied
    ! that has previously been used for something else:
    nullify(mesh%faces)
    nullify(mesh%columns)
    nullify(mesh%element_columns)

    allocate(mesh%colourings(NUM_COLOURINGS))
    do i = 1, NUM_COLOURINGS
       nullify(mesh%colourings(i)%sets)
    end do
    allocate(mesh%ndglno(elements*shape%loc))

#ifdef _OPENMP
    ! Use first touch policy.
    !$OMP PARALLEL DO SCHEDULE(STATIC)
    do i=1, mesh%elements
       do j=1, shape%loc
          mesh%ndglno((i-1)*shape%loc+j)=0
       end do
    end do
    !$OMP END PARALLEL DO
#endif

#ifdef HAVE_MEMORY_STATS
    call register_allocation("mesh_type", "integer", elements*shape%loc,&
         & name=mesh%name)
#endif

    allocate(mesh%adj_lists)
    mesh%wrapped=.false.
    nullify(mesh%region_ids)
    nullify(mesh%subdomain_mesh)
    nullify(mesh%refcount) ! Hack for gfortran component initialisation
    !                         bug.
    mesh%periodic=.false.
    
    call addref(mesh)

  end subroutine allocate_mesh

! IAKOVOS commented out
!  subroutine allocate_scalar_field(field, mesh, name, field_type, py_func, py_positions)

! IAKOVOS commented out
!  subroutine allocate_vector_field(field, dim, mesh, name, field_type)

! IAKOVOS commented out
!  subroutine allocate_tensor_field(field, mesh, name, field_type, dim)
  
  subroutine deallocate_subdomain_mesh(mesh)
    type(mesh_type) :: mesh

    if (.not.associated(mesh%subdomain_mesh)) return

    deallocate(mesh%subdomain_mesh%element_list)
    deallocate(mesh%subdomain_mesh%node_list)

    deallocate(mesh%subdomain_mesh)

  end subroutine deallocate_subdomain_mesh

! IAKOVOS commented out
!  subroutine deallocate_mesh_faces(mesh)

! IAKOVOS commented out
!  subroutine deallocate_mesh(mesh)

! IAKOVOS commented out
!  recursive subroutine deallocate_scalar_field(field)
    
! IAKOVOS commented out
!  subroutine remove_boundary_conditions_scalar(field)
 
! IAKOVOS commented out
!  recursive subroutine deallocate_vector_field(field)
 
! IAKOVOS commented out
!  subroutine remove_boundary_conditions_vector(field)
  
! IAKOVOS commented out
!  subroutine deallocate_tensor_field(field)
  
! IAKOVOS commented out
!  subroutine allocate_scalar_boundary_condition(bc, mesh, surface_element_list, &
!    name, type)
    
! IAKOVOS commented out
!  subroutine allocate_vector_boundary_condition(bc, mesh, surface_element_list, &
    
! IAKOVOS commented out
!  subroutine deallocate_scalar_boundary_condition(bc)
  
! IAKOVOS commented out
!  subroutine deallocate_vector_boundary_condition(bc)
    
  !---------------------------------------------------------------------
  ! routines for wrapping meshes and fields around provided arrays
  !---------------------------------------------------------------------
  
  function wrap_mesh(ndglno, shape, name) result (mesh)
    !!< Return a mesh wrapped around the information provided.
    type(mesh_type) :: mesh

    integer, dimension(:), target, intent(in) :: ndglno
    type(element_type), target, intent(in) :: shape
    character(len=*), intent(in) :: name

    mesh%ndglno=>ndglno
    mesh%shape=shape
    call incref(shape)
    nullify(mesh%faces)

    mesh%name=name
    
    mesh%elements=size(ndglno)/shape%loc

    allocate(mesh%adj_lists)
    mesh%wrapped=.true.
    mesh%nodes=maxval(ndglno)
    nullify(mesh%refcount) ! Hack for gfortran component initialisation
    !                         bug.
    mesh%periodic = .false. ! can only really assume that this is false as
                            ! we have no other information
    call addref(mesh)

  end function wrap_mesh

! IAKOVOS commented out
!  function wrap_scalar_field(mesh, val, name, val_stride) result (field)

  ! IAKOVOS commented out
!  function wrap_tensor_field(mesh, val, name) result (field)

! IAKOVOS commented out
!  function make_mesh (model, shape, continuity, name) &
!       result (mesh)

  ! IAKOVOS commented out
!  subroutine add_faces(mesh, model, sndgln, sngi, boundary_ids, &
!    periodic_face_map, element_owner, incomplete_surface_mesh, stat)

  ! IAKOVOS commented out
!  subroutine add_faces_face_list(mesh, sndgln, boundary_ids, &
!    element_owner, incomplete_surface_mesh)

! IAKOVOS commented out
!  subroutine register_internal_surface_element(mesh, sele, ele, neighbour_ele)

! IAKOVOS commented out
!  subroutine register_external_surface_element(mesh, sele, ele, snodes)

  ! IAKOVOS commented out
!  subroutine add_faces_face_list_periodic_from_non_periodic_model( &
!     mesh, model, periodic_face_map)
    
  ! IAKOVOS commented out
!  subroutine add_faces_face_list_non_periodic_from_periodic_model( &
!     mesh, model, periodic_face_map, stat)
    
! IAKOVOS commented out
!  subroutine fix_periodic_face_orientation(nonperiodic, periodic, periodic_face_map)

! IAKOVOS commented out
!  subroutine create_surface_mesh(surface_mesh, surface_nodes, &
!    mesh, surface_elements, name)
    
  logical function SetContains(a, b)
  !!< Auxillary function that returns true if b contains a
  integer, dimension(:), intent(in):: a, b
  
    integer i
    
    SetContains=.false.
    do i=1, size(a)
      if (.not. any(b==a(i))) return
    end do
    SetContains=.true.

  end function SetContains

! IAKOVOS commented out
!  function make_mesh_periodic(positions,physical_boundary_ids,aliased_boundary_ids,periodic_mapping_python,name, &
!       periodic_face_map) result (positions_out)

  ! IAKOVOS commented out
!  function make_fake_mesh_linearnonconforming(model, name) result (mesh)

  ! IAKOVOS commented out
!  function make_submesh (model, name) &
!       result (mesh)

  ! IAKOVOS commented out
!  function extract_elements(positions, elements) result(subpos)
  
  ! IAKOVOS commented out
!  subroutine add_lists_mesh(mesh, nnlist, nelist, eelist)
  
! IAKOVOS commented out
!  subroutine add_lists_scalar(field, nnlist, nelist, eelist)
  
! IAKOVOS commented out
!  subroutine add_lists_vector(field, nnlist, nelist, eelist)
  
! IAKOVOS commented out
!  subroutine add_lists_tensor(field, nnlist, nelist, eelist)

  ! IAKOVOS commented out
!  subroutine extract_lists_mesh(mesh, nnlist, nelist, eelist)
  
  ! IAKOVOS commented out
!  subroutine extract_lists_scalar(field, nnlist, nelist, eelist)
 
  ! IAKOVOS commented out
!  subroutine extract_lists_vector(field, nnlist, nelist, eelist)

  ! IAKOVOS commented out
!  subroutine extract_lists_tensor(field, nnlist, nelist, eelist)
  
  ! IAKOVOS commented out
!  subroutine add_nnlist_mesh(mesh)
  
! IAKOVOS commented out
!  subroutine add_nnlist_scalar(field)
 
! IAKOVOS commented out
!  subroutine add_nnlist_vector(field)
 
! IAKOVOS commented out
!  subroutine add_nnlist_tensor(field)
  
  ! IAKOVOS commented out
!  function extract_nnlist_mesh(mesh) result(nnlist)
 
  ! IAKOVOS commented out 
!  function extract_nnlist_scalar(field) result(nnlist)
 
  ! IAKOVOS commented out
!  function extract_nnlist_vector(field) result(nnlist)
  
  ! IAKOVOS commented out
!  function extract_nnlist_tensor(field) result(nnlist)

  ! IAKOVOS commented out
!  subroutine add_nelist_mesh(mesh)
 
  ! IAKOVOS commented out
!  subroutine add_nelist_scalar(field)
 
  ! IAKOVOS commented out
!  subroutine add_nelist_vector(field)
  
  ! IAKOVOS commented out
!  subroutine add_nelist_tensor(field)
  
  ! IAKOVOS commented out
!  function extract_nelist_mesh(mesh) result(nelist)
 
  ! IAKOVOS commented out
!  function extract_nelist_scalar(field) result(nelist)
 
  ! IAKOVOS commented out
!  function extract_nelist_vector(field) result(nelist)
  
  ! IAKOVOS commented out
!  function extract_nelist_tensor(field) result(nelist)

  ! IAKOVOS commented out
!  subroutine add_eelist_mesh(mesh)

  ! IAKOVOS commented out
!  subroutine add_eelist_scalar(field)
  
  ! IAKOVOS commented out
!  subroutine add_eelist_vector(field)
  
  ! IAKOVOS commented out
!  subroutine add_eelist_tensor(field)

  function extract_eelist_mesh(mesh) result(eelist)
    !!< Extract the element-element list (generating if necessary) from the
    !!< adjacency cache for the supplied mesh
  
    type(mesh_type), intent(in) :: mesh
    
    type(csr_sparsity), pointer :: eelist
    
    call add_eelist(mesh)
    eelist => mesh%adj_lists%eelist
    assert(has_references(eelist))
    
  end function extract_eelist_mesh
 
  ! IAKOVOS commented out
!  function extract_eelist_scalar(field) result(eelist)
 
  function extract_eelist_vector(field) result(eelist)
    !!< Extract the element-element list (generating if necessary) from the
    !!< adjacency cache for the supplied field
    
    type(vector_field), intent(in) :: field
    
    type(csr_sparsity), pointer :: eelist
    
    eelist => extract_eelist(field%mesh)
    
  end function extract_eelist_vector
 
  ! IAKOVOS commented out
!  function extract_eelist_tensor(field) result(eelist)
  
  subroutine remove_lists_mesh(mesh)
    !!< Remove the adjecency lists from the adjacency cache for the supplied
    !!< mesh
  
    type(mesh_type), intent(inout) :: mesh
    
    call remove_nnlist(mesh)
    call remove_nelist(mesh)
    call remove_eelist(mesh)
  
  end subroutine remove_lists_mesh
  
  ! IAKOVOS commented out
!  subroutine remove_nnlist_mesh(mesh)

  ! IAKOVOS commented out
!  subroutine remove_nelist_mesh(mesh)

  ! IAKOVOS commented out
!  subroutine remove_eelist_mesh(mesh)
  
  subroutine zero_scalar(field)
    !!< Set all entries in the field provided to 0.0
    type(scalar_field), intent(inout) :: field
#ifdef _OPENMP
    integer :: i
#endif
    
    assert(field%field_type/=FIELD_TYPE_PYTHON)
    
#ifdef _OPENMP
    ! Use first touch policy.
    !$OMP PARALLEL DO SCHEDULE(STATIC)
    do i=1, size(field%val)
       field%val(i)=0.0
    end do
    !$OMP END PARALLEL DO
#else
    field%val=0.0
#endif

  end subroutine zero_scalar

  subroutine zero_vector(field)
    !!< Set all entries in the field provided to 0.0
    type(vector_field), intent(inout) :: field

#ifdef _OPENMP
    integer :: i
#endif

    assert(field%field_type/=FIELD_TYPE_PYTHON)
    
#ifdef _OPENMP
    ! Use first touch policy.
    !$OMP PARALLEL DO SCHEDULE(STATIC)
    do i=1, size(field%val, 2)
       field%val(:,i)=0.0
    end do
    !$OMP END PARALLEL DO
#else
       field%val=0.0
#endif

  end subroutine zero_vector

  subroutine zero_vector_dim(field, dim)
    !!< Set all entries in dimension dim of the field provided to 0.0
    type(vector_field), intent(inout) :: field
    integer, intent(in) :: dim

#ifdef _OPENMP
    integer :: j
#endif

    assert(field%field_type/=FIELD_TYPE_PYTHON)

#ifdef _OPENMP
       ! Use first touch policy.
       !$OMP PARALLEL DO SCHEDULE(STATIC)
       do j=1, size(field%val, 2)
          field%val(dim,j)=0.0
       end do
       !$OMP END PARALLEL DO
#else
       field%val(dim,:)=0.0
#endif

  end subroutine zero_vector_dim

  subroutine zero_tensor(field)
    !!< Set all entries in the field provided to 0.0
    type(tensor_field), intent(inout) :: field

#ifdef _OPENMP
    integer :: j
#endif

    assert(field%field_type/=FIELD_TYPE_PYTHON)
    
#ifdef _OPENMP
    ! Use first touch policy.
    !$OMP PARALLEL DO SCHEDULE(STATIC)
    do j=1, size(field%val, 3)
       field%val(:,:,j)=0.0
    end do
    !$OMP END PARALLEL DO
#else
    field%val=0.0
#endif

  end subroutine zero_tensor  

  subroutine zero_tensor_dim_dim(field, dim1, dim2)
    !!< Set all entries in the component indicated of field to 0.0
    type(tensor_field), intent(inout) :: field
    integer, intent(in) :: dim1, dim2

#ifdef _OPENMP
    integer :: j
#endif

    assert(field%field_type/=FIELD_TYPE_PYTHON)

#ifdef _OPENMP
    ! Use first touch policy.
    !$OMP PARALLEL DO SCHEDULE(STATIC)
    do j=1, size(field%val, 3)
       field%val(dim1,dim2,j)=0.0
    end do
    !$OMP END PARALLEL DO
#else
    field%val(dim1,dim2,:)=0.0
#endif
    
  end subroutine zero_tensor_dim_dim

  subroutine zero_scalar_field_nodes(field, node_numbers)
    !!< Zeroes the scalar field at the specified node_numbers
    !!< Does not work for constant fields
    type(scalar_field), intent(inout) :: field
    integer, dimension(:), intent(in) :: node_numbers

    assert(field%field_type==FIELD_TYPE_NORMAL)
    
    field%val(node_numbers) = 0.0
    
  end subroutine zero_scalar_field_nodes
  
  subroutine zero_vector_field_nodes(field, node_numbers)
    !!< Zeroes the vector field at the specified nodes
    !!< Does not work for constant fields
    type(vector_field), intent(inout) :: field
    integer, dimension(:), intent(in) :: node_numbers
    integer :: i

    assert(field%field_type==FIELD_TYPE_NORMAL)
    
    do i=1,field%dim
      field%val(i,node_numbers) = 0.0
    end do
    
  end subroutine zero_vector_field_nodes

  subroutine zero_tensor_field_nodes(field, node_numbers)
    !!< Zeroes the tensor field at the specified nodes
    !!< Does not work for constant fields
    type(tensor_field), intent(inout) :: field
    integer, dimension(:), intent(in) :: node_numbers

    assert(field%field_type==FIELD_TYPE_NORMAL)

    field%val(:, :, node_numbers) = 0.0
    
  end subroutine zero_tensor_field_nodes

#include "Reference_count_mesh_type.F90"
#include "Reference_count_scalar_field.F90"
#include "Reference_count_vector_field.F90"
#include "Reference_count_tensor_field.F90"

end module libsupermesh_fields_allocates
