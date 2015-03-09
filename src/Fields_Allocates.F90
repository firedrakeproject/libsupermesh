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
use libsupermesh_shape_functions, only: make_element_shape
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
use libsupermesh_ieee_arithmetic
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
  public :: create_surface_mesh
! IAKOVOS commented out
!  public :: extract_scalar_field, wrap_mesh, wrap_scalar_field, &
!    & wrap_tensor_field
  public :: wrap_mesh, make_mesh
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
     module procedure allocate_scalar_field, allocate_vector_field, &
           & allocate_mesh
  end interface

  interface deallocate
! IAKOVOS commented out
!     module procedure deallocate_mesh, deallocate_scalar_field,&
!          & deallocate_vector_field, deallocate_tensor_field, &
!          & deallocate_scalar_boundary_condition, &
!          & deallocate_vector_boundary_condition
     module procedure deallocate_mesh
  end interface

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

  subroutine allocate_scalar_field(field, mesh, name, field_type, py_func, py_positions)
    type(scalar_field), intent(out) :: field
    type(mesh_type), intent(in), target :: mesh
    character(len=*), intent(in),optional :: name
    integer, intent(in), optional :: field_type

    character(len=*), intent(in), optional ::  py_func
    type(vector_field), intent(in), optional, target :: py_positions

    integer :: lfield_type
    integer :: stat
    integer :: toloc, fromloc

    if (present(field_type)) then
      lfield_type = field_type
    else
      lfield_type = FIELD_TYPE_NORMAL
    end if

    field%mesh=mesh
    call incref(mesh)
    
    if (present(name)) then
       field%name=name
    else
       field%name=empty_name
    end if

    field%field_type = lfield_type
    select case(lfield_type)
    case(FIELD_TYPE_NORMAL)
      allocate(field%val(node_count(mesh)))
      field%py_dim = mesh_dim(mesh)
      field%py_positions_shape => mesh%shape

#ifdef HAVE_MEMORY_STATS
      call register_allocation("scalar_field", "real", node_count(mesh), &
        name=name)
#endif
    case(FIELD_TYPE_CONSTANT)
      allocate(field%val(1))
      field%py_dim = mesh_dim(mesh)
      field%py_positions_shape => mesh%shape

#ifdef HAVE_MEMORY_STATS
      call register_allocation("scalar_field", "real", 1, name=name)
#endif
    case(FIELD_TYPE_DEFERRED)
      allocate(field%val(0))
      field%py_dim = mesh_dim(mesh)
      field%py_positions_shape => mesh%shape
    case(FIELD_TYPE_PYTHON)
      if (present(py_func)) then
        field%py_func = py_func
      else
        if (stat /= 0) then
!          FLAbort("Field specified as FIELD_TYPE_PYTHON, but no func passed!")		! ToDo
        end if
      end if

      if (.not. present(py_positions)) then
!        FLAbort("Field specified as FIELD_TYPE_PYTHON but no positions field passed!")	! ToDo
      end if
      field%py_positions => py_positions
      field%py_dim = py_positions%dim
      field%py_positions_shape => py_positions%mesh%shape
      call incref(field%py_positions_shape)
      call incref(field%py_positions)

      if (associated(py_positions%mesh%refcount, mesh%refcount)) then
        field%py_positions_same_mesh = .true.
      else
        field%py_positions_same_mesh = .false.
        allocate(field%py_locweight(mesh%shape%loc, py_positions%mesh%shape%loc))
        do toloc=1,size(field%py_locweight,1)
          do fromloc=1,size(field%py_locweight,2)
            field%py_locweight(toloc,fromloc)=eval_shape(py_positions%mesh%shape, fromloc, &
              local_coords(toloc, mesh%shape))
          end do
        end do
      end if

      call add_nelist(field%mesh)
    end select

    field%wrapped=.false.
    field%aliased=.false.
    field%option_path=empty_path
    allocate(field%bc)
    nullify(field%refcount) ! Hacks for gfortran component initialisation
    !                         bug.
    call addref(field)

    call zero(field)
    
  end subroutine allocate_scalar_field

  subroutine allocate_vector_field(field, dim, mesh, name, field_type)
    type(vector_field), intent(out) :: field
    integer, intent(in) :: dim
    type(mesh_type), intent(in), target :: mesh
    character(len=*), intent(in), optional :: name
    integer, intent(in), optional :: field_type
    integer :: n_count
    integer :: lfield_type

    if (present(field_type)) then
      lfield_type = field_type
    else
      lfield_type = FIELD_TYPE_NORMAL
    end if
    
    field%dim=dim
    field%option_path=empty_path

    field%mesh=mesh
    call incref(mesh)
    
    if (present(name)) then
       field%name=name
    else
       field%name=empty_name
    end if

    field%field_type = lfield_type
    select case(lfield_type)
    case(FIELD_TYPE_NORMAL)
      n_count = node_count(mesh)
      allocate(field%val(dim,n_count))
#ifdef HAVE_MEMORY_STATS
      call register_allocation("vector_field", "real", n_count*dim, &
        name=name)
#endif
    case(FIELD_TYPE_CONSTANT)
      allocate(field%val(dim,1))
#ifdef HAVE_MEMORY_STATS
      call register_allocation("vector_field", "real", dim, name=name)
#endif
    case(FIELD_TYPE_DEFERRED)
      allocate(field%val(0,0))
    end select

    field%wrapped = .false.
    field%aliased = .false.
    allocate(field%bc)
    nullify(field%refcount) ! Hack for gfortran component initialisation
    !                         bug.    
    
    allocate(field%picker)
    
    call addref(field)

    call zero(field)

  end subroutine allocate_vector_field

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

  subroutine deallocate_mesh(mesh)
    !!< Deallocate the components of mesh. Shape functions are not
    !!< deallocated here.
    type(mesh_type), intent(inout) :: mesh
    integer :: i
    call decref(mesh)
    if (has_references(mesh)) then
       ! There are still references to this mesh so we don't deallocate.
       return
    end if
    call deallocate(mesh%shape)
    
    if (.not.mesh%wrapped) then
#ifdef HAVE_MEMORY_STATS
       call register_deallocation("mesh_type", "integer", &
            size(mesh%ndglno), name=mesh%name)
#endif
       deallocate(mesh%ndglno)
    end if

    if(associated(mesh%region_ids)) then
       deallocate(mesh%region_ids)
    end if
    
    assert(associated(mesh%adj_lists))
    call remove_lists(mesh)
    deallocate(mesh%adj_lists)
    nullify(mesh%adj_lists)
    
    if(associated(mesh%halos)) then
       call deallocate(mesh%halos)
       deallocate(mesh%halos)
    end if

    if(associated(mesh%element_halos)) then
       call deallocate(mesh%element_halos)
       deallocate(mesh%element_halos)
    end if

    call deallocate_faces(mesh)

    if(associated(mesh%subdomain_mesh)) then
       call deallocate_subdomain_mesh(mesh)
    end if
    
    if(associated(mesh%columns)) then
      deallocate(mesh%columns)
    end if
    
    if(associated(mesh%element_columns)) then
      deallocate(mesh%element_columns)
    end if

    if(associated(mesh%colourings)) then
       do i = 1, NUM_COLOURINGS
          if(associated(mesh%colourings(i)%sets)) then
             call deallocate(mesh%colourings(i)%sets)
             deallocate(mesh%colourings(i)%sets)
          end if
       end do
       deallocate(mesh%colourings)
    end if
  end subroutine deallocate_mesh

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

  function make_mesh (model, shape, continuity, name) &
       result (mesh)
    !!< Produce a mesh based on an old mesh but with a different shape and/or continuity.
    type(mesh_type) :: mesh

    type(mesh_type), intent(in) :: model
    type(element_type), target, intent(in), optional :: shape
    integer, intent(in), optional :: continuity
    character(len=*), intent(in), optional :: name
    
    integer, dimension(:), allocatable :: ndglno
    real, dimension(:), pointer :: val
    integer :: i, input_nodes, n_faces
#ifdef _OPENMP
    integer :: j
#endif

    if (present(continuity)) then
       mesh%continuity=continuity
    else
       mesh%continuity=model%continuity
    end if

    allocate(mesh%adj_lists)
    mesh%elements=model%elements
    mesh%periodic=model%periodic
    mesh%wrapped=.false.

    if (present(shape)) then
       mesh%shape=shape
    else
       mesh%shape=model%shape
    end if
    call incref(mesh%shape)

    ! You can't have a CG degree 0 mesh!
    if(mesh%shape%degree==0.and.mesh%continuity>=0.and.mesh%shape&
         &%numbering%type/=ELEMENT_TRACE) then
!      FLExit("For a P0 mesh, the 'mesh_continuity' must be Discontinuous.")		! ToDo
    end if

    if (present(name)) then
       mesh%name=name
    else
       mesh%name=empty_name
    end if

    if (associated(model%region_ids)) then
       allocate(mesh%region_ids(size(model%region_ids)))
       mesh%region_ids=model%region_ids
    end if

    if (mesh%continuity>=0) then
       ! Make a continuous field.
       if (model%continuity<0) then
!          FLExit("Unable to derive a continuous mesh from a discontinuous mesh")	! ToDo
       end if

       allocate(ndglno(mesh%shape%numbering%vertices*model%elements), &
            mesh%ndglno(mesh%shape%loc*model%elements))

#ifdef _OPENMP
          ! Use first touch policy.
          !$OMP PARALLEL DO SCHEDULE(STATIC)
          do i=1, mesh%elements
             do j=1, mesh%shape%loc
                mesh%ndglno((i-1)*mesh%shape%loc+j)=0
             end do
          end do
          !$OMP END PARALLEL DO
#endif

#ifdef HAVE_MEMORY_STATS
       call register_allocation("mesh_type", "integer", &
            size(mesh%ndglno), name=name)
#endif

       if(model%shape%degree==1 .or. ele_count(model) == 0) then
          ndglno=model%ndglno
          input_nodes = node_count(model)
       else
          ndglno=mesh_connectivity(model)
          input_nodes = maxval(ndglno)
       end if
       
       if (associated(model%halos)) then
          assert(element_halo_count(model) > 0)
          allocate(mesh%halos(size(model%halos)))

          call make_global_numbering &
               (mesh%nodes, mesh%ndglno, input_nodes, mesh%elements, &
               ndglno, mesh%shape, model%halos, model%element_halos(1), &
               mesh%halos)

          allocate(mesh%element_halos(size(model%element_halos)))
          do i=1,size(mesh%element_halos)
             mesh%element_halos(i)=model%element_halos(i)
             call incref(mesh%element_halos(i))
          end do

          do i=1,size(mesh%halos)
             call reorder_halo_from_element_halo(mesh%halos(i), mesh&
                  &%element_halos(1), mesh)
          end do

       else
          
          call make_global_numbering &
               (mesh%nodes, mesh%ndglno, max(maxval(ndglno), 0), mesh%elements, &
               ndglno, mesh%shape)
       end if

    else
       !trace fields have continuity -1 but aren't like DG
       if(mesh%shape%numbering%type/=ELEMENT_TRACE) then
          ! Make a discontinuous field.
          allocate(mesh%ndglno(mesh%shape%loc*model%elements))

#ifdef _OPENMP
          ! Use first touch policy.
          !$OMP PARALLEL DO SCHEDULE(STATIC)
          do i=1, mesh%elements
             do j=1, mesh%shape%loc
                mesh%ndglno((i-1)*mesh%shape%loc+j)=0
             end do
          end do
          !$OMP END PARALLEL DO
#endif

#ifdef HAVE_MEMORY_STATS
          call register_allocation("mesh_type", "integer", &
               size(mesh%ndglno), name=name)
#endif
          if (associated(model%halos)) then
             assert(associated(model%element_halos))
             allocate(mesh%halos(size(model%halos)))
             
             
             call make_global_numbering_DG(mesh%nodes, mesh%ndglno, &
                  mesh%elements, mesh%shape, model%element_halos, &
                  mesh%halos)
             
             allocate(mesh%element_halos(size(model%element_halos)))
             do i=1,size(mesh%element_halos)
                mesh%element_halos(i)=model%element_halos(i)
                call incref(mesh%element_halos(i))
             end do
             
          else
             
             call make_global_numbering_DG(mesh%nodes, mesh%ndglno, &
                  mesh%elements, mesh%shape)
             
          end if
       end if
    end if

    nullify(mesh%refcount) ! Hack for gfortran component initialisation
    !                         bug.
    
    ! Transfer the eelist from model to mesh
    assert(associated(model%adj_lists))
    if(associated(model%adj_lists%eelist)) then
      ewrite(2, *) "Transferring element-element list to mesh " // trim(mesh%name)
      allocate(mesh%adj_lists%eelist)
      mesh%adj_lists%eelist = model%adj_lists%eelist
      call incref(mesh%adj_lists%eelist)
    end if
    
    if(has_faces(model)) then
      call add_faces(mesh, model)
    end if

    if (mesh%shape%numbering%type==ELEMENT_TRACE) then
       select case(mesh%shape%numbering%family)
       case(FAMILY_SIMPLEX)          
          n_faces = mesh%shape%dim + 1
       case(FAMILY_CUBE)
          n_faces = 2*mesh%shape%dim
       case default
!          FLExit('Element family not supported for trace elements')		! ToDo
       end select
       allocate(mesh%ndglno(mesh%elements*n_faces*mesh%faces%shape%loc))
       call make_global_numbering_trace(mesh)
       call create_surface_mesh(mesh%faces%surface_mesh, &
            mesh%faces%surface_node_list, mesh, name='Surface'//trim(mesh%name))
#ifdef HAVE_MEMORY_STATS
       call register_allocation("mesh_type", "integer", &
            size(mesh%faces%surface_node_list), name='Surface'//trim(mesh%name))
#endif
    end if
    call addref(mesh)

  end function make_mesh

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

  subroutine create_surface_mesh(surface_mesh, surface_nodes, &
    mesh, surface_elements, name)
  !! Creates a surface mesh consisting of the surface elements
  !! specified by surface_element_list  
  type(mesh_type), intent(out):: surface_mesh
  !! Returns a pointer to a list containing the global node number
  !! of the nodes on this surface mesh, can be used for surface node
  !! to global node numbering conversion.
  integer, dimension(:), pointer:: surface_nodes
  !! mesh to take surface from (should have %faces component)
  type(mesh_type), intent(in):: mesh
  !! which surface elements to select
  !! (if not provided all surface elements are included)
  integer, dimension(:), optional, target,intent(in):: surface_elements
  !! name for the new surface_mesh
  character(len=*), intent(in):: name
  
    integer, dimension(:), pointer:: lsurface_elements
    integer, dimension(:), pointer:: suf_ndglno
    integer, dimension(:), allocatable:: nod2sufnod
    integer, dimension(mesh%faces%shape%loc):: glnodes
    integer i, j, sele, sufnod, snloc
    
    snloc=mesh%faces%shape%loc
    
    if (present(surface_elements)) then
       lsurface_elements => surface_elements
    else
       allocate(lsurface_elements(1:surface_element_count(mesh)))
       lsurface_elements=(/ (i, i=1, size(lsurface_elements)) /)
    end if
      
    allocate(nod2sufnod(1:node_count(mesh)))
    nod2sufnod=0
    
    ! mark surface nodes with nod2sufnod(nod)==1
    sufnod=0
    do i=1, size(lsurface_elements)
      sele=lsurface_elements(i)
      glnodes=face_global_nodes(mesh, sele)
      do j=1, snloc
        if (nod2sufnod(glnodes(j))==0) then
          sufnod=sufnod+1
          nod2sufnod(glnodes(j))=1
        end if
      end do
    end do
      
    call allocate(surface_mesh, nodes=sufnod, &
      elements=size(lsurface_elements), shape=mesh%faces%shape, &
      name=name)
      
    surface_mesh%periodic=mesh%periodic
    surface_mesh%continuity=mesh%continuity

    allocate(surface_nodes(1:sufnod))
    
    ! create numbering in the same order as full nodal numbering:
    sufnod=0
    do i=1, size(nod2sufnod)
      if (nod2sufnod(i)==1) then
        sufnod=sufnod+1
        ! global node to surface node numbering
        nod2sufnod(i)=sufnod
        ! and the reverse
        surface_nodes(sufnod)=i
        
      end if
    end do
      
    ! map global node numbering to surface node numbering
    suf_ndglno => surface_mesh%ndglno
    do i=1, size(lsurface_elements)
      sele=lsurface_elements(i)
      suf_ndglno( (i-1)*snloc+1:i*snloc )=nod2sufnod(face_global_nodes(mesh, sele))
    end do
        
    deallocate(nod2sufnod)
    
    if (.not. present(surface_elements)) then
       deallocate(lsurface_elements)
    end if
  
  end subroutine create_surface_mesh
    
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
