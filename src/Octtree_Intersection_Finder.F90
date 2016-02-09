#include "fdebug.h"

#define TREE_DIM 3
#define TREE_NCHILDREN 8

module libsupermesh_octtree_intersection_finder

  use libsupermesh_fldebug, only : flabort_pinpoint
  use libsupermesh_fldebug_parameters, only : debug_log_unit
  use libsupermesh_integer_set, only : integer_set, insert
  use libsupermesh_intersections, only : intersections, deallocate, &
    & intersections_to_csr_sparsity
  use libsupermesh_quadtree_intersection_finder, only : max_nelist_degree

  implicit none

  private
  
  public :: octtree_node, octtree_element, max_nelist_degree, &
    & octtree_intersection_finder, build_octtree, query_octtree, deallocate
  
  ! Octtree node
  type octtree_node
    ! Number of stored elements. If this is <= max_n then this is a leaf node.
    integer :: n = 0
    ! Maximum number of stored elements
    integer :: max_n  
    ! Node bounding box
    real, dimension(2, TREE_DIM) :: bbox
    ! If this is a leaf node, elements stored by this node
    type(octtree_element), dimension(:), pointer :: elements
    ! Otherwise, children of this node
    type(octtree_node), dimension(:), pointer :: children => null()
  end type octtree_node
  
  ! Octtree stored elements
  type octtree_element
    ! Element bounding box
    real, dimension(2, TREE_DIM) :: bbox
    ! Element index
    integer :: ele
  end type octtree_element
  
  interface octtree_intersection_finder
    module procedure octtree_intersection_finder_intersections, &
      & octtree_intersection_finder_csr_sparsity
  end interface octtree_intersection_finder
  
  interface query_octtree
    module procedure query_octtree_internal, query_octtree_integer_set
  end interface query_octtree
  
  interface deallocate
    module procedure deallocate_node
  end interface deallocate
    
contains

  subroutine octtree_intersection_finder_intersections(positions_a, enlist_a, positions_b, enlist_b, map_ab, max_size)
    ! dim x nnodes_a
    real, dimension(:, :), intent(in) :: positions_a
    ! loc_a x nelements_a
    integer, dimension(:, :), intent(in) :: enlist_a
    ! dim x nnodes_b
    real, dimension(:, :), intent(in) :: positions_b
    ! loc_b x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    ! nelements_a
    type(intersections), dimension(:), intent(out) :: map_ab
    integer, optional, intent(in) :: max_size
    
    integer :: ele_a, nelements_a, nelements_b
    type(octtree_node) :: root_node
    
    integer :: neles_b
    integer, dimension(:), allocatable :: eles_b
    logical, dimension(:), allocatable :: seen_ele_b
    
    if(size(positions_a, 1) /= TREE_DIM) then
      FLAbort("Invalid dimension")
    end if
    nelements_a = size(enlist_a, 2)
    nelements_b = size(enlist_b, 2)
    if(nelements_a == 0 .or. nelements_b == 0) then
      do ele_a = 1, nelements_a
        allocate(map_ab(ele_a)%v(0))
        map_ab(ele_a)%n = 0
      end do
      return
    end if
    
    root_node = build_octtree(positions_b, enlist_b, max_size = max_size)
    
    allocate(eles_b(nelements_b), seen_ele_b(nelements_b))   
    neles_b = 0
    seen_ele_b = .false. 
    do ele_a = 1, nelements_a
      seen_ele_b(eles_b(:neles_b)) = .false.
      neles_b = 0
      call query_octtree(root_node, bbox(positions_a(:, enlist_a(:, ele_a))), eles_b, neles_b, seen_ele_b)
      map_ab(ele_a)%n = neles_b
      allocate(map_ab(ele_a)%v(neles_b))
      map_ab(ele_a)%v = eles_b(:neles_b)
    end do
    deallocate(eles_b, seen_ele_b)
    
    call deallocate(root_node)
    
  end subroutine octtree_intersection_finder_intersections
  
  subroutine octtree_intersection_finder_csr_sparsity(positions_a, enlist_a, positions_b, enlist_b, map_ab_indices, map_ab_indptr, max_size)
    ! dim x nnodes_a
    real, dimension(:, :), intent(in) :: positions_a
    ! loc_a x nelements_a
    integer, dimension(:, :), intent(in) :: enlist_a
    ! dim x nnodes_b
    real, dimension(:, :), intent(in) :: positions_b
    ! loc_b x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    ! Compressed Sparse Row (CSR) sparsity pattern, as described in:
    !   http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.sparse.csr_matrix.html
    integer, dimension(:), allocatable, intent(out) :: map_ab_indices
    ! nelements_a + 1
    integer, dimension(:), intent(out) :: map_ab_indptr
    integer, optional, intent(in) :: max_size
    
    type(intersections), dimension(:), allocatable :: map_ab
    
    allocate(map_ab(size(enlist_a, 2)))
    call octtree_intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab, max_size = max_size)
    call intersections_to_csr_sparsity(map_ab, map_ab_indices, map_ab_indptr)
    call deallocate(map_ab)
    deallocate(map_ab)
    
  end subroutine octtree_intersection_finder_csr_sparsity

  function build_octtree(positions, enlist, max_size) result(root_node)
    ! TREE_DIM x nnodes
    real, dimension(:, :), intent(in) :: positions
    ! loc x nelements
    integer, dimension(:, :), intent(in) :: enlist
    integer, optional, intent(in) :: max_size
    
    type(octtree_node) :: root_node
    
    integer, parameter :: default_max_size = 256 
    
    integer :: ele, nelements, node, nnodes
    
    if(size(positions, 1) /= TREE_DIM) then
      FLAbort("Invalid dimension")
    end if
    nnodes = size(positions, 2)
    if(nnodes == 0) then
      FLAbort("No nodes")
    end if
    if(size(enlist, 1) == 0) then
      FLAbort("No element local nodes")
    end if
    nelements = size(enlist, 2)
    
    if(present(max_size)) then
      ! A maximum number of stored elements has been provided. It is the
      ! caller's responsibility to ensure that this is large enough.
      root_node%max_n = max_size
    else
      root_node%max_n = max(max_nelist_degree(nnodes, enlist), default_max_size)
    end if
    allocate(root_node%elements(root_node%max_n))
    
    ! Compute the root node bounding box
    root_node%bbox(1, :) = positions(:, 1)
    root_node%bbox(2, :) = positions(:, 1)
    do node = 2, nnodes
      root_node%bbox(1, 1) = min(root_node%bbox(1, 1), positions(1, node))
      root_node%bbox(1, 2) = min(root_node%bbox(1, 2), positions(2, node))
      root_node%bbox(1, 3) = min(root_node%bbox(1, 3), positions(3, node))
      root_node%bbox(2, 1) = max(root_node%bbox(2, 1), positions(1, node))
      root_node%bbox(2, 2) = max(root_node%bbox(2, 2), positions(2, node))
      root_node%bbox(2, 3) = max(root_node%bbox(2, 3), positions(3, node))
    end do    
    
    ! Add all elements
    do ele = 1, nelements
      call add_element(root_node, ele, bbox(positions(:, enlist(:, ele))))
    end do
    
  end function build_octtree
  
  pure recursive subroutine deallocate_node(node)
    type(octtree_node), intent(inout) :: node
    
    integer :: i
    
    if(node%n <= node%max_n) then
      deallocate(node%elements)
    else
      do i = 1, TREE_NCHILDREN
        call deallocate(node%children(i))
      end do
      deallocate(node%children)
      nullify(node%children)
    end if
    node%n = 0
    
  end subroutine deallocate_node
  
  pure recursive subroutine add_element(node, ele, bbox)
    type(octtree_node), intent(inout) :: node
    integer, intent(in) :: ele
    real, dimension(2, TREE_DIM), intent(in) :: bbox
  
    integer :: i, j
    real, dimension(TREE_DIM) :: mid_point
  
    if(.not. bboxes_intersect(node%bbox, bbox)) return
    
    if(node%n <= node%max_n - 1) then
      ! Add new element to the leaf node
      node%n = node%n + 1
      node%elements(node%n)%bbox = bbox
      node%elements(node%n)%ele = ele
    else if(node%n == node%max_n) then
      ! Current leaf node is full. Split and add elements to the new child leaf
      ! nodes.
    
      ! Create new child leaf nodes
      allocate(node%children(TREE_NCHILDREN))
      mid_point = 0.5 * (node%bbox(1, :) + node%bbox(2, :))
      node%children(1)%bbox(1, :) = (/node%bbox(1, 1), node%bbox(1, 2), node%bbox(1, 3)/)
      node%children(1)%bbox(2, :) = (/mid_point(1),    mid_point(2),    mid_point(3)/)      
      node%children(2)%bbox(1, :) = (/node%bbox(1, 1), mid_point(2),    node%bbox(1, 3)/)
      node%children(2)%bbox(2, :) = (/mid_point(1),    node%bbox(2, 2), mid_point(3)/)
      node%children(3)%bbox(1, :) = (/mid_point(1),    mid_point(2),    node%bbox(1, 3)/)
      node%children(3)%bbox(2, :) = (/node%bbox(2, 1), node%bbox(2, 2), mid_point(3)/)
      node%children(4)%bbox(1, :) = (/mid_point(1),    node%bbox(1, 2), node%bbox(1, 3)/)
      node%children(4)%bbox(2, :) = (/node%bbox(2, 1), mid_point(2),    mid_point(3)/)
      node%children(5)%bbox(1, :) = (/node%bbox(1, 1), node%bbox(1, 2), mid_point(3)/)
      node%children(5)%bbox(2, :) = (/mid_point(1),    mid_point(2),    node%bbox(2, 3)/)      
      node%children(6)%bbox(1, :) = (/node%bbox(1, 1), mid_point(2),    mid_point(3)/)
      node%children(6)%bbox(2, :) = (/mid_point(1),    node%bbox(2, 2), node%bbox(2, 3)/)
      node%children(7)%bbox(1, :) = (/mid_point(1),    mid_point(2),    mid_point(3)/)
      node%children(7)%bbox(2, :) = (/node%bbox(2, 1), node%bbox(2, 2), node%bbox(2, 3)/)
      node%children(8)%bbox(1, :) = (/mid_point(1),    node%bbox(1, 2), mid_point(3)/)
      node%children(8)%bbox(2, :) = (/node%bbox(2, 1), mid_point(2),    node%bbox(2, 3)/)
      do i = 1, TREE_NCHILDREN
        node%children(i)%max_n = node%max_n
        allocate(node%children(i)%elements(node%max_n))
      end do
      
      ! Add elements stored in the parent to the new child leaf nodes
      do i = 1, node%max_n
        do j = 1, TREE_NCHILDREN
          call add_element(node%children(j), node%elements(i)%ele, node%elements(i)%bbox)
        end do
      end do
      ! Add new element to the new child leaf nodes
      do j = 1, TREE_NCHILDREN
        call add_element(node%children(j), ele, bbox)
      end do
      
      ! Mark the parent as not a leaf node
      node%n = huge(0)
      deallocate(node%elements)
      nullify(node%elements)
    else
      ! Not a leaf node. Add to child nodes.
      do j = 1, TREE_NCHILDREN
        call add_element(node%children(j), ele, bbox)
      end do
    end if
  
  end subroutine add_element
  
  pure recursive subroutine query_octtree_internal(node, bbox_a, eles_b, neles_b, seen_ele_b)
    type(octtree_node), intent(in) :: node
    real, dimension(2, TREE_DIM), intent(in) :: bbox_a
    integer, dimension(:), intent(inout) :: eles_b
    integer, intent(inout) :: neles_b
    logical, dimension(:), intent(inout) :: seen_ele_b
    
    integer :: ele_b, i
    
    if(node%n == 0) return
    if(.not. bboxes_intersect(node%bbox, bbox_a)) return
    
    if(node%n <= node%max_n) then
      ! A leaf node. Return stored elements whose bounding boxes intersect with
      ! the query bounding box.
      do i = 1, node%n
        ele_b = node%elements(i)%ele
        if(.not. seen_ele_b(ele_b)) then
          if(bboxes_intersect(node%elements(i)%bbox, bbox_a)) then
            neles_b = neles_b + 1
            eles_b(neles_b) = ele_b
            seen_ele_b(ele_b) = .true.
          end if
        end if
      end do
    else
      ! Not a leaf node. Query child nodes.
      do i = 1, TREE_NCHILDREN
        call query_octtree(node%children(i), bbox_a, eles_b, neles_b, seen_ele_b)
      end do
    end if
    
  end subroutine query_octtree_internal
  
  recursive subroutine query_octtree_integer_set(node, bbox_a, eles_b)
    type(octtree_node), intent(in) :: node
    real, dimension(2, TREE_DIM), intent(in) :: bbox_a
    type(integer_set), intent(inout) :: eles_b
    
    integer :: ele_b, i
    
    if(node%n == 0) return
    if(.not. bboxes_intersect(node%bbox, bbox_a)) return
    
    if(node%n <= node%max_n) then
      ! A leaf node. Return stored elements whose bounding boxes intersect with
      ! the query bounding box.
      do i = 1, node%n
        ele_b = node%elements(i)%ele
        if(bboxes_intersect(node%elements(i)%bbox, bbox_a)) then
          call insert(eles_b, ele_b)
        end if
      end do
    else
      ! Not a leaf node. Query child nodes.
      do i = 1, TREE_NCHILDREN
        call query_octtree(node%children(i), bbox_a, eles_b)
      end do
    end if
    
  end subroutine query_octtree_integer_set

  pure function bbox(coords)
    ! TREE_DIM x loc
    real, dimension(:, :), intent(in) :: coords

    real, dimension(2, TREE_DIM) :: bbox

    integer :: i

    bbox(1, :) = coords(:, 1)
    bbox(2, :) = coords(:, 1)
    do i = 2, size(coords, 2)
      bbox(1, 1) = min(bbox(1, 1), coords(1, i))
      bbox(1, 2) = min(bbox(1, 2), coords(2, i))
      bbox(1, 3) = min(bbox(1, 3), coords(3, i))
      bbox(2, 1) = max(bbox(2, 1), coords(1, i))
      bbox(2, 2) = max(bbox(2, 2), coords(2, i))
      bbox(2, 3) = max(bbox(2, 3), coords(3, i))
    end do

  end function bbox

  pure function bboxes_intersect(bbox_1, bbox_2) result(intersect)
    real, dimension(2, TREE_DIM), intent(in) :: bbox_1
    real, dimension(2, TREE_DIM), intent(in) :: bbox_2

    logical :: intersect

    intersect = bbox_2(2, 1) > bbox_1(1, 1) .and. bbox_2(1, 1) < bbox_1(2, 1) &
        & .and. bbox_2(2, 2) > bbox_1(1, 2) .and. bbox_2(1, 2) < bbox_1(2, 2) &
        & .and. bbox_2(2, 3) > bbox_1(1, 3) .and. bbox_2(1, 3) < bbox_1(2, 3)

  end function bboxes_intersect

end module libsupermesh_octtree_intersection_finder
