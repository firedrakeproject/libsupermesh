#include "fdebug.h"

module libsupermesh_fields_dummy

  use libsupermesh_fldebug

  implicit none

  private

  public :: mesh_type, vector_field, eelist_type, allocate, deallocate, &
    & mesh_eelist, sloc, extract_eelist, ele_count, node_count, ele_val, &
    & ele_loc, node_val, row_m_ptr, ele_nodes, set, set_ele_nodes, &
    & simplex_volume, triangle_area, tetrahedron_volume

  type mesh_type
    integer :: dim
    integer :: nnodes
    integer :: nelements
    integer :: loc
    integer :: sloc
    integer, dimension(:), pointer :: ndglno
    integer :: continuity
    type(eelist_ptr), pointer :: eelist
    
    integer, pointer :: refcount
  end type mesh_type

  type vector_field
    type(mesh_type) :: mesh
    integer :: dim
    real, dimension(:, :), pointer :: val
    
    integer, pointer :: refcount
  end type vector_field

  type eelist_type
    integer, dimension(:, :), pointer :: v
    integer, dimension(:), pointer :: n
  end type eelist_type

  type eelist_ptr
    type(eelist_type), pointer :: ptr
  end type eelist_ptr

  interface allocate
    module procedure allocate_mesh, allocate_vector_field
  end interface allocate

  interface deallocate
    module procedure deallocate_mesh, deallocate_vector_field, deallocate_eelist
  end interface deallocate

  interface extract_eelist
    module procedure extract_eelist_mesh, extract_eelist_vector_field
  end interface extract_eelist

  interface ele_count
    module procedure ele_count_vector_field
  end interface ele_count

  interface node_count
    module procedure node_count_vector_field
  end interface node_count

  interface ele_val
    module procedure ele_val_vector_field
  end interface ele_val

  interface ele_loc
    module procedure ele_loc_vector_field
  end interface ele_loc

  interface node_val
    module procedure node_val_vector_field
  end interface node_val

  interface row_m_ptr
    module procedure row_m_ptr_eelist
  end interface row_m_ptr

  interface ele_nodes
    module procedure ele_nodes_vector_field
  end interface ele_nodes

  interface set
    module procedure set_vector_field_nodes, set_vector_field_node
  end interface set

  interface set_ele_nodes
    module procedure set_ele_nodes_mesh
  end interface set_ele_nodes

contains

  pure subroutine allocate_mesh(mesh, dim, nnodes, nelements, loc)
    type(mesh_type), intent(out) :: mesh
    integer, intent(in) :: dim
    integer, intent(in) :: nnodes
    integer, intent(in) :: nelements
    integer, intent(in) :: loc

    mesh%dim = dim
    mesh%nnodes = nnodes
    mesh%nelements = nelements
    mesh%loc = loc
    mesh%sloc = sloc(dim, loc)
    allocate(mesh%ndglno(nelements * loc))
    allocate(mesh%eelist)
    nullify(mesh%eelist%ptr)
    mesh%continuity = 0
    
    allocate(mesh%refcount)
    mesh%refcount = 1
    
  end subroutine allocate_mesh

  pure subroutine deallocate_mesh(mesh)
    type(mesh_type), intent(inout) :: mesh

    mesh%refcount = mesh%refcount - 1
    if(mesh%refcount == 0) then

      deallocate(mesh%ndglno)
      if(associated(mesh%eelist%ptr)) then
        call deallocate(mesh%eelist%ptr)
        deallocate(mesh%eelist%ptr)
      end if
      deallocate(mesh%eelist)

      deallocate(mesh%refcount)
    end if

  end subroutine deallocate_mesh

  subroutine allocate_vector_field(field, dim, mesh)
    type(vector_field), intent(out) :: field
    integer, intent(in) :: dim
    type(mesh_type), intent(in) :: mesh

    field%dim = dim
    allocate(field%val(dim, mesh%nnodes))
    field%mesh = mesh

    mesh%refcount = mesh%refcount + 1
    allocate(field%refcount)
    field%refcount = 1
    
  end subroutine allocate_vector_field

  pure subroutine deallocate_vector_field(field)
    type(vector_field), intent(inout) :: field

    field%refcount = field%refcount - 1
    if(field%refcount == 0) then
    
      call deallocate(field%mesh)
      deallocate(field%val)
      
      deallocate(field%refcount)      
    end if

  end subroutine deallocate_vector_field

  subroutine mesh_eelist(nnodes, enlist, sloc, eelist)
    integer, intent(in) :: nnodes
    ! loc x nelements
    integer, dimension(:, :), intent(in) :: enlist
    integer, intent(in) :: sloc
    type(eelist_type), intent(out) :: eelist

    integer :: ele, i, j, k, k_max, k_min, loc, nelements, nneigh_ele

    integer, dimension(:), allocatable :: nelist_n
    ! Compressed Sparse Row (CSR) sparsity pattern, as described in:
    !   http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.sparse.csr_matrix.html
    integer, dimension(:), allocatable :: nelist_indices, nelist_indptr

    integer, dimension(:), allocatable :: seen, seen_eles
    integer :: nseen_eles

    loc = size(enlist, 1)
    nelements = size(enlist, 2)

    ! Construct the node-element list. Based on NODELE in Fluidity
    ! femtools/Adjacency_Lists.F90 (see e.g. Fluidity version 4.1.11). Code
    ! first added 05/07/2015.
    allocate(nelist_n(nnodes), nelist_indptr(nnodes + 1))
    nelist_n = 0
    do i = 1, nelements
      do j = 1, loc
        nelist_n(enlist(j, i)) = nelist_n(enlist(j, i)) + 1
      end do
    end do
    nelist_indptr(1) = 1
    do i = 1, nnodes
      nelist_indptr(i + 1) = nelist_indptr(i) + nelist_n(i)
    end do
    allocate(nelist_indices(nelist_indptr(nnodes + 1) - 1))
    nelist_n = 0
    do i = 1, nelements
      do j = 1, loc
        nelist_indices(nelist_indptr(enlist(j, i)) + nelist_n(enlist(j, i))) = i
        nelist_n(enlist(j, i)) = nelist_n(enlist(j, i)) + 1
      end do
    end do
    deallocate(nelist_n)

    allocate(seen_eles(nelements), seen(nelements))
    seen = 0
    nneigh_ele = loc  ! Linear simplex assumption here (works for more general
                      ! cases, but less efficient)
    allocate(eelist%v(nneigh_ele, nelements), eelist%n(nelements))
    eelist%n = 0
    do i = 1, nelements
      nseen_eles = 0
      loc_loop: do j = 1, loc
        k_min = nelist_indptr(enlist(j, i))
        k_max = nelist_indptr(enlist(j, i) + 1) - 1
        do k = k_min, k_max
          ele = nelist_indices(k)
          if(ele /= i) then
            seen(ele) = seen(ele) + 1
            if(seen(ele) == 1) then
              nseen_eles = nseen_eles + 1
              seen_eles(nseen_eles) = ele
            end if
            if(seen(ele) == sloc) then
              eelist%n(i) = eelist%n(i) + 1
#ifndef NDEBUG
              if(eelist%n(i) > nneigh_ele) then
                FLAbort("Invalid connectivity")
              end if
#endif
              eelist%v(eelist%n(i), i) = ele
#ifdef NDEBUG
              if(eelist%n(i) == nneigh_ele) exit loc_loop
#else
            else if(seen(ele) > sloc) then
              FLAbort("Invalid connectivity")
#endif
            end if
          end if
        end do
      end do loc_loop
      do j = 1, nseen_eles
        seen(seen_eles(j)) = 0
      end do
    end do
    deallocate(nelist_indices, nelist_indptr, seen_eles, seen)

  end subroutine mesh_eelist

  pure function sloc(dim, loc)
    integer, intent(in) :: dim
    integer, intent(in) :: loc

    integer :: sloc

    if(loc == dim + 1) then
      ! Simplex
      sloc = loc - 1
    else if(loc == 2 ** dim) then
      ! Cubic
      sloc = loc / 2
    else
      sloc = -1
    end if

  end function sloc

  pure subroutine deallocate_eelist(eelist)
    type(eelist_type), intent(inout) :: eelist

    deallocate(eelist%v, eelist%n)

  end subroutine deallocate_eelist

  function extract_eelist_mesh(mesh) result(eelist)
    type(mesh_type), intent(inout) :: mesh

    type(eelist_type), pointer :: eelist

    if(.not. associated(mesh%eelist%ptr)) then
      allocate(mesh%eelist%ptr)
      call mesh_eelist(mesh%nnodes, reshape(mesh%ndglno, (/mesh%loc, mesh%nelements/)), mesh%sloc, mesh%eelist%ptr)
    end if
    eelist => mesh%eelist%ptr

  end function extract_eelist_mesh

  function extract_eelist_vector_field(field) result(eelist)
    type(vector_field), intent(inout) :: field

    type(eelist_type), pointer :: eelist

    eelist => extract_eelist(field%mesh)

  end function extract_eelist_vector_field

  pure function ele_count_vector_field(field) result(ele_count)
    type(vector_field), intent(in) :: field

    integer :: ele_count

    ele_count = field%mesh%nelements

  end function ele_count_vector_field

  pure function node_count_vector_field(field) result(node_count)
    type(vector_field), intent(in) :: field

    integer :: node_count

    node_count = field%mesh%nnodes

  end function node_count_vector_field

  pure function ele_val_vector_field(field, ele) result(ele_val)
    type(vector_field), intent(in) :: field
    integer, intent(in) :: ele

    real, dimension(field%dim, field%mesh%loc) :: ele_val

    integer :: i, i_0

    i_0 = (ele - 1) * field%mesh%loc
    do i = 1, field%mesh%loc
      ele_val(:, i) = field%val(:, field%mesh%ndglno(i_0 + i))
    end do

  end function ele_val_vector_field

  pure function ele_loc_vector_field(field, ele) result(ele_loc)
    type(vector_field), intent(in) :: field
    integer, intent(in) :: ele

    integer :: ele_loc

    ele_loc = field%mesh%loc

  end function ele_loc_vector_field

  pure function node_val_vector_field(field, node) result(node_val)
    type(vector_field), intent(in) :: field
    integer, intent(in) :: node

    real, dimension(field%dim) :: node_val

    node_val = field%val(:, node)

  end function node_val_vector_field

  function row_m_ptr_eelist(eelist, i) result(row_m_ptr)
    type(eelist_type), intent(in) :: eelist
    integer, intent(in) :: i

    integer, dimension(:), pointer :: row_m_ptr

    row_m_ptr => eelist%v(:eelist%n(i), i)

  end function row_m_ptr_eelist

  function ele_nodes_vector_field(field, ele) result(nodes)
    type(vector_field), intent(in) :: field
    integer, intent(in) :: ele

    integer, dimension(:), pointer :: nodes

    nodes => field%mesh%ndglno((ele - 1) * field%mesh%loc + 1:ele * field%mesh%loc)

  end function ele_nodes_vector_field

  pure subroutine set_vector_field_nodes(field, nodes, val)
    type(vector_field), intent(inout) :: field
    integer, dimension(:), intent(in) :: nodes
    ! dim x loc
    real, dimension(:, :), intent(in) :: val

    field%val(:, nodes) = val

  end subroutine set_vector_field_nodes

  pure subroutine set_vector_field_node(field, node, val)
    type(vector_field), intent(inout) :: field
    integer, intent(in) :: node
    ! dim
    real, dimension(:), intent(in) :: val

    field%val(:, node) = val

  end subroutine set_vector_field_node

  pure subroutine set_ele_nodes_mesh(mesh, ele, nodes)
    type(mesh_type), intent(inout) :: mesh
    integer, intent(in) :: ele
    ! loc
    integer, dimension(:), intent(in) :: nodes

    mesh%ndglno((ele - 1) * mesh%loc + 1:ele * mesh%loc) = nodes

  end subroutine set_ele_nodes_mesh

  pure function simplex_volume(cell_coords) result(volume)
    ! dim x loc
    real, dimension(:, :), intent(in) :: cell_coords

    real :: volume

    select case(size(cell_coords, 1))
      case(1)
        volume = abs(cell_coords(1, 2) - cell_coords(1, 1))
      case(2)
        volume = triangle_area(cell_coords)
      case(3)
        volume = tetrahedron_volume(cell_coords)
      case default
        volume = -huge(0.0)
    end select
  
  end function simplex_volume

  pure function triangle_area(cell_coords) result(area)
    real, dimension(2, 3), intent(in) :: cell_coords

    real :: area

    real, dimension(2) :: e1, e2

    e1 = cell_coords(:, 2) - cell_coords(:, 1)
    e2 = cell_coords(:, 3) - cell_coords(:, 1)

    area = 0.5 * abs(e1(1) * e2(2) - e1(2) * e2(1))

  end function triangle_area

  pure function tetrahedron_volume(cell_coords) result(volume)
    real, dimension(3, 4), intent(in) :: cell_coords

    real :: volume

    real, dimension(3) :: e1, e2, e3

    e1 = cell_coords(:, 2) - cell_coords(:, 1)
    e2 = cell_coords(:, 3) - cell_coords(:, 1)
    e3 = cell_coords(:, 4) - cell_coords(:, 1)

    volume = (1.0 / 6.0) * abs(e1(1) * (e2(2) * e3(3) - e2(3) * e3(2)) &
                           & + e1(2) * (e2(3) * e3(1) - e2(1) * e3(3)) &
                           & + e1(3) * (e2(1) * e3(2) - e2(2) * e3(1)))
    
  end function tetrahedron_volume

end module libsupermesh_fields_dummy
