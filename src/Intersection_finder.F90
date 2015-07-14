!#define COUNT_INTERSECTION_TESTS
#include "fdebug.h"

module libsupermesh_intersection_finder_module

 use libsupermesh_fields_dummy
 use libsupermesh_integer_hash_table
 use libsupermesh_integer_set
 use libsupermesh_linked_lists

implicit none

interface crtree_intersection_finder_set_input
  subroutine libsupermesh_cintersection_finder_set_input(positions, enlist, ndim, loc, nnodes, nelements) bind(c)
    use iso_c_binding
    implicit none
    integer(kind = c_int), intent(in) :: ndim, loc, nnodes, nelements
    real(kind = c_double), intent(in), dimension(nnodes * ndim) :: positions
    integer(kind = c_int), intent(in), dimension(nelements * loc) :: enlist
  end subroutine libsupermesh_cintersection_finder_set_input
end interface crtree_intersection_finder_set_input

interface crtree_intersection_finder_find
  subroutine libsupermesh_cintersection_finder_find(positions, ndim, loc) bind(c)
    use iso_c_binding
    implicit none
    integer(kind = c_int), intent(in) :: ndim, loc
    real(kind = c_double), dimension(ndim * loc) :: positions
  end subroutine libsupermesh_cintersection_finder_find
end interface crtree_intersection_finder_find

interface rtree_intersection_finder_query_output
  subroutine libsupermesh_cintersection_finder_query_output(nelems) bind(c)
    use iso_c_binding
    implicit none
    integer(kind = c_int), intent(out) :: nelems
  end subroutine libsupermesh_cintersection_finder_query_output
end interface rtree_intersection_finder_query_output

interface rtree_intersection_finder_get_output
  subroutine libsupermesh_cintersection_finder_get_output(id, nelem) bind(c)
    use iso_c_binding
    implicit none
    integer(kind = c_int), intent(out) :: id
    integer(kind = c_int), intent(in) :: nelem
  end subroutine libsupermesh_cintersection_finder_get_output
end interface rtree_intersection_finder_get_output

interface rtree_intersection_finder_reset
  subroutine libsupermesh_cintersection_finder_reset(ntests) bind(c)
    use iso_c_binding
    implicit none
    integer(kind = c_int), intent(out) :: ntests
  end subroutine libsupermesh_cintersection_finder_reset
end interface rtree_intersection_finder_reset

private

public :: rtree_intersection_finder_set_input, rtree_intersection_finder_find, &
  & rtree_intersection_finder_query_output, &
  & rtree_intersection_finder_get_output, rtree_intersection_finder_reset

public :: bbox_predicate, intersection_tests, &
  & reset_intersection_tests_counter, intersection_finder, &
  & advancing_front_intersection_finder_seeds, advancing_front_intersection_finder, &
  & rtree_intersection_finder, brute_force_intersection_finder
  
#ifdef COUNT_INTERSECTION_TESTS
integer, save :: intersection_tests_count = 0
#endif

contains

  function bbox(pos) result(box)
    ! dim x loc
    real, dimension(:, :) :: pos
    real, dimension(size(pos, 1), 2) :: box
    integer :: dim, i, loc, j

    dim = size(pos, 1)
    loc = size(pos, 2)
    do i=1,dim
      box(i, 1) = pos(i, 1)
      box(i, 2) = pos(i, 1)
      do j=2,loc
        box(i, 1) = min(pos(i, j), box(i, 1))
        box(i, 2) = max(pos(i, j), box(i, 2))
      end do
    end do
    
  end function bbox

  function bbox_predicate(bboxA, bboxB) result(intersects)
    real, dimension(:, :), intent(in) :: bboxA, bboxB
    logical :: intersects
    integer :: dim, i

#ifdef COUNT_INTERSECTION_TESTS
    intersection_tests_count = intersection_tests_count + 1
#endif

    intersects = .false.
    dim = size(bboxA, 1)

    do i=1,dim
      if (bboxA(i, 1) > bboxB(i, 2)) return
      if (bboxA(i, 2) < bboxB(i, 1)) return
    end do
    intersects = .true.

  end function bbox_predicate
  
  function intersection_tests() result(tests)
    !!< Return the number of intersection tests
    
    integer :: tests
    
#ifdef COUNT_INTERSECTION_TESTS
    tests = intersection_tests_count
#else
    FLAbort("Counting of intersection tests is not available")
    ! To keep the compiler quiet
    tests = 0
#endif

  end function intersection_tests
  
  subroutine reset_intersection_tests_counter()
    !!< Reset the intersection tests counter
    
#ifdef COUNT_INTERSECTION_TESTS
    intersection_tests_count = 0
#else
    FLAbort("Counting of intersection tests is not available")
#endif

  end subroutine reset_intersection_tests_counter
      
  function intersection_finder(positionsA, ndglnoA, positionsB, ndglnoB, seed) result(map_AB)
    ! The positions and meshes of A and B
    real, dimension(:, :), intent(in) :: positionsA
    integer, dimension(:, :), intent(in) :: ndglnoA
    real, dimension(:, :), intent(in) :: positionsB
    integer, dimension(:, :), intent(in) :: ndglnoB
    ! for each element in A, the intersecting elements in B
    type(ilist), dimension(size(ndglnoA, 2)) :: map_AB
    integer, optional, intent(in) :: seed
    !!< A simple wrapper to select an intersection finder

    integer :: i, dimA, dimB, nodesA, nodesB, verticesA, verticesB
    type(mesh_type) :: mesh
    type(vector_field), target :: lpositionsA
    
    type(ilist) :: seeds
    type(inode), pointer :: node
    type(ilist), dimension(:), allocatable :: sub_map_AB

    ewrite(1, *) "In intersection_finder"

    ! We cannot assume connectedness, so we may have to run the
    ! advancing front more than once (once per connected sub-domain)

    dimA = size(positionsA, 1)
    nodesA = size(positionsA, 2)
    dimB = size(positionsB, 1)
    nodesB = size(positionsB, 2)
    if(.not. dimA == dimB) then
      FLAbort("Incompatible input meshes")
    end if
    verticesA = size(ndglnoA, 1)
    verticesB = size(ndglnoB, 1)
    
    call allocate(mesh, dimA, nodesA, size(ndglnoA, 2), size(ndglnoA, 1))
    do i = 1, size(ndglnoA, 2)
      mesh%ndglno((i - 1) * verticesA + 1:i * verticesA) = ndglnoA(:, i)
    end do
    call allocate(lpositionsA, dimA, mesh)
    call deallocate(mesh)
    lpositionsA%val = positionsA
    
    seeds = advancing_front_intersection_finder_seeds(lpositionsA)
    
    call deallocate(lpositionsA)

    allocate(sub_map_AB(size(map_AB)))
    node => seeds%firstnode
    do while(associated(node))
      sub_map_AB = advancing_front_intersection_finder( &
        & positionsA, ndglnoA, &
        & positionsB, ndglnoB, &
        & seed = node%value)
      do i = 1, size(sub_map_AB)
        if(sub_map_AB(i)%length > 0) then
          assert(map_AB(i)%length == 0)
          map_AB(i) = sub_map_AB(i)
        end if
      end do
    
      node => node%next
    end do
    
    deallocate(sub_map_AB)
    call deallocate(seeds)

    ewrite(1, *) "Exiting intersection_finder"
  
  end function intersection_finder
  
  function connected(positions)
    !!< Return whether the supplied coordinate field is connected. Uses a simple
    !!< element advancing front.
    
    type(vector_field), intent(inout) :: positions
    
    logical :: connected
    
    integer :: ele, i
    type(eelist_type), pointer :: eelist
    logical, dimension(:), allocatable :: tested
    integer, dimension(:), pointer :: neigh
    type(ilist) :: next
    
    eelist => extract_eelist(positions)
    
    allocate(tested(ele_count(positions)))
    tested = .false.
    
    assert(ele_count(positions) > 0)
    ele = 1
    tested(ele) = .true.
    neigh => row_m_ptr(eelist, ele)
    do i = 1, size(neigh)
      if(neigh(i) <= 0) cycle
    
      call insert(next, neigh(i))
    end do
    
    do while(next%length > 0)
      ele = pop(next)
      if(tested(ele)) cycle
      
      tested(ele) = .true.
      neigh => row_m_ptr(eelist, ele)
      do i = 1, size(neigh)
        if(neigh(i) <= 0) cycle
        if(tested(neigh(i))) cycle
        ! Should check if neigh(i) is already in the list
        
        call insert(next, neigh(i))
      end do
    end do
    
    ! The mesh is connected iff we see all elements in the advancing front
    connected = all(tested)
    
    deallocate(tested)
    
  end function connected
  
  function advancing_front_intersection_finder_seeds(positions) result(seeds)
    !!< Return a list of seeds for the advancing front intersection finder - one
    !!< seed per connected sub-domain.

    type(vector_field), intent(inout) :: positions

    type(ilist) :: seeds
    integer :: ele, first_ele, i
    type(eelist_type), pointer :: eelist
    logical, dimension(:), allocatable :: tested
    integer, dimension(:), pointer :: neigh
    type(ilist) :: next

    eelist => extract_eelist(positions)
    
    allocate(tested(ele_count(positions)))
    tested = .false.
    
    first_ele = 1
    do while(first_ele /= 0)
      ele = first_ele
      assert(ele > 0)
      assert(ele <= ele_count(positions))
      assert(.not. tested(ele))
      
      call insert(seeds, ele)
      
      tested(ele) = .true.
      neigh => row_m_ptr(eelist, ele)
      do i = 1, size(neigh)
        if(neigh(i) <= 0) cycle
        if(tested(neigh(i))) cycle
      
        call insert(next, neigh(i))
      end do
      
      do while(next%length > 0)
        ele = pop(next)
        if(tested(ele)) cycle
        
        tested(ele) = .true.
        neigh => row_m_ptr(eelist, ele)
        do i = 1, size(neigh)
          if(neigh(i) <= 0) cycle
          if(tested(neigh(i))) cycle
          ! Should check if neigh(i) is already in the list
          
          call insert(next, neigh(i))
        end do
      end do
      
      first_ele = next_false_loc(first_ele + 1, tested)
    end do
    assert(all(tested))

    deallocate(tested)
    
  contains
  
    pure function next_false_loc(start_index, logical_vector) result(loc)
      integer, intent(in) :: start_index
      logical, dimension(:), intent(in) :: logical_vector

      integer :: loc, i

      do i = start_index, size(logical_vector)
        if(.not. logical_vector(i)) then
          loc = i
          return
        end if
      end do
      
      loc = 0
      
    end function next_false_loc
    
  end function advancing_front_intersection_finder_seeds

  function advancing_front_intersection_finder(positionsA, ndglnoA, positionsB, ndglnoB, seed) result(map_AB)
    real, dimension(:, :), intent(in) :: positionsA
    integer, dimension(:, :), intent(in) :: ndglnoA
    real, dimension(:, :), intent(in) :: positionsB
    integer, dimension(:, :), intent(in) :: ndglnoB
    integer, optional, intent(in) :: seed

    type(ilist), dimension(size(ndglnoA, 2)) :: map_AB

    integer :: dimA, dimB, nodesA, nodesB, verticesA, verticesB
    type(mesh_type) :: mesh
    type(vector_field), target :: lpositionsA, lpositionsB

    ! processed_neighbour maps an element to a neighbour that has already been processed (i.e. its clue)
    type(integer_hash_table) :: processed_neighbour
    ! we also need to keep a set of the elements we've seen: this is different to
    ! the elements that have map_AB(ele)%length > 0 in the case where the domain
    ! is not simply connected!
    type(integer_set) :: seen_elements

    integer :: ele_A
    type(mesh_type), pointer :: mesh_A, mesh_B
    integer :: i, neighbour
    real, dimension(size(ndglnoB, 2), size(positionsB, 1), 2) :: bboxes_B
    integer, dimension(:), pointer :: neigh_A
    type(eelist_type), pointer :: eelist_A, eelist_B

    type(ilist) :: clues

    ewrite(1, *) "In advancing_front_intersection_finder"

    dimA = size(positionsA, 1)
    nodesA = size(positionsA, 2)
    dimB = size(positionsB, 1)
    nodesB = size(positionsB, 2)
    if(.not. dimA == dimB) then
      FLAbort("Incompatible input meshes")
    end if
    verticesA = size(ndglnoA, 1)
    verticesB = size(ndglnoB, 1)

    call allocate(mesh, dimA, nodesA, size(ndglnoA, 2), size(ndglnoA, 1))
    do i = 1, size(ndglnoA, 2)
      mesh%ndglno((i - 1) * verticesA + 1:i * verticesA) = ndglnoA(:, i)
    end do
    call allocate(lpositionsA, dimA, mesh)
    call deallocate(mesh)
    lpositionsA%val = positionsA

    call allocate(mesh, dimB, nodesB, size(ndglnoB, 2), size(ndglnoB, 1))
    do i = 1, size(ndglnoB, 2)
      mesh%ndglno((i - 1) * verticesB + 1:i * verticesB) = ndglnoB(:, i)
    end do
    call allocate(lpositionsB, dimB, mesh)
    call deallocate(mesh)
    lpositionsB%val = positionsB

    mesh_A => lpositionsA%mesh
    mesh_B => lpositionsB%mesh

    eelist_A => extract_eelist(mesh_A)
    eelist_B => extract_eelist(mesh_B)

    call compute_bboxes(lpositionsB, bboxes_B)

    if(present(seed)) then
      assert(seed > 0)
      assert(seed <= ele_count(lpositionsA))
      ele_A = seed
    else
      ele_A = 1
    end if
    map_AB(ele_A) = brute_force_search(ele_val(lpositionsA, ele_A), bboxes_B)

    call allocate(processed_neighbour)
    call allocate(seen_elements)

    neigh_A => row_m_ptr(eelist_A, ele_A)
    do i=1,size(neigh_A)
      neighbour = neigh_A(i)
      if (neighbour <= 0) cycle
      call insert(processed_neighbour, neighbour, ele_A)
    end do
    call insert(seen_elements, ele_A)

    do while (key_count(processed_neighbour) > 0)
      call fetch_pair(processed_neighbour, 1, ele_A, neighbour)
      ! try to keep our memory footprint low
      call remove(processed_neighbour, ele_A)
      call insert(seen_elements, ele_A)

      assert(map_AB(ele_A)%length == 0) ! we haven't seen it yet
      clues = clueful_search(ele_val(lpositionsA, ele_A), map_AB(neighbour), &
                           & bboxes_B, ele_A, neighbour)
      map_AB(ele_A) = advance_front(ele_val(lpositionsA, ele_A), lpositionsB, clues, bboxes_B, eelist_B)
      call deallocate(clues)

      ! Now that ele_A has been computed, make its clues available to anyone who needs them
      neigh_A => row_m_ptr(eelist_A, ele_A)
      do i=1,size(neigh_A)
        neighbour = neigh_A(i)
        if (neighbour <= 0) cycle
        if (has_value(seen_elements, neighbour)) then
          ! We've already seen it
          cycle
        end if
        call insert(processed_neighbour, neighbour, ele_A)
      end do
    end do

    assert(key_count(processed_neighbour) == 0)
    call deallocate(processed_neighbour)
    call deallocate(seen_elements)
    call deallocate(lpositionsA)
    call deallocate(lpositionsB)

    ewrite(1, *) "Exiting advancing_front_intersection_finder"

    contains
      function advance_front(posA, positionsB, clues, bboxes_B, eelist_B) result(map)
        real, dimension(:, :), intent(in) :: posA
        type(vector_field), intent(in), target :: positionsB
        type(ilist), intent(inout) :: clues
        real, dimension(:, :, :), intent(in) :: bboxes_B
        type(eelist_type), intent(in) :: eelist_B

        type(ilist) :: map
        integer, dimension(:), pointer :: neigh_B
        integer :: i, possible, neighbour, j
        logical :: intersects
        type(mesh_type), pointer :: mesh_B
        real, dimension(size(posA, 1), 2) :: bboxA
        integer :: ele_B
        type(integer_set) :: in_list
        type(integer_hash_table) :: possibles_tbl
        integer :: possible_size

        bboxA = bbox(posA)
        call allocate(in_list)
        call allocate(possibles_tbl)
        possible_size = 0

        mesh_B => positionsB%mesh

        do while (clues%length /= 0)
          ele_B = pop(clues)
          if (.not. has_value(in_list, ele_B)) then
            call insert(map, ele_B)
            call insert(in_list, ele_B)
          end if

          ! Append all the neighbours of ele_B to possibles. 
          neigh_B => row_m_ptr(eelist_B, ele_B)
          do i=1,size(neigh_B)
            neighbour = neigh_B(i)
            if (neighbour <= 0) cycle
            if (.not. has_value(in_list, neighbour)) then
              possible_size = possible_size + 1
              call insert(possibles_tbl, possible_size, neighbour)
              call insert(in_list, neighbour)
            end if
          end do
        end do

        ! while len(possibles) != 0:
          ! If predicate(ele_A, ele_B) is false: remove it from possibles and add it to rejects.
          ! If true: add it to map and add all its neighbours not in map or rejects to map.

        j = 1
        do while (j <= possible_size)
          possible = fetch(possibles_tbl, j)
          intersects = bbox_predicate(bboxA, bboxes_B(possible, :, :))
          if (intersects) then
            call insert(map, possible)
            neigh_B => row_m_ptr(eelist_B, possible)
            do i=1,size(neigh_B)
              neighbour = neigh_B(i)
              if (neighbour <= 0) cycle
              if (.not. has_value(in_list, neighbour)) then
                possible_size = possible_size + 1
                call insert(possibles_tbl, possible_size, neighbour)
                call insert(in_list, neighbour)
              end if
            end do
          end if
          j = j + 1
        end do

        call deallocate(in_list)
        call deallocate(possibles_tbl)

        possible_size = 0
      end function advance_front
  end function advancing_front_intersection_finder
  
  function brute_force_intersection_finder(positions_a, positions_b) result(map_ab)
    !!< As advancing_front_intersection_finder, but uses a brute force
    !!< algorithm. For testing *only*. For practical applications, use the
    !!< linear algorithm.
  
    ! The positions and meshes of A and B
    type(vector_field), intent(in), target :: positions_a, positions_b
    ! for each element in A, the intersecting elements in B
    type(ilist), dimension(ele_count(positions_a)) :: map_ab
    
    integer :: i, j
    
    ewrite(1, *) "In brute_force_intersection_finder"
    
    do i = 1, ele_count(positions_a)
      do j = 1, ele_count(positions_b)
        if(bbox_predicate(bbox(ele_val(positions_a, i)), bbox(ele_val(positions_b, j)))) then
          call insert(map_ab(i), j)
        end if
      end do
    end do
    
    ewrite(1, *) "Exiting brute_force_intersection_finder"
    
  end function brute_force_intersection_finder
  
!  subroutine rtree_intersection_finder_set_input(positions, ndim, nnodes, elementCount, loc, enlist)
  subroutine rtree_intersection_finder_set_input(old_positions)
    type(vector_field), intent(in) :: old_positions
!    real, intent(in), dimension(nnodes * ndim) :: positions
!    integer, intent(in) :: ndim, nnodes, elementCount, loc
!    integer, intent(in), dimension(elementCount * loc) :: enlist
    real, dimension(node_count(old_positions) * old_positions%dim) :: tmp_positions
    integer :: node, dim

    dim = old_positions%dim

    ! Ugh. We have to copy the memory because old_positions
    ! stores it as 2 or 3 separate vectors
    do node=1,node_count(old_positions)
      tmp_positions((node-1)*dim+1:node*dim) = node_val(old_positions, node)
    end do

    call crtree_intersection_finder_set_input(tmp_positions, old_positions%mesh%ndglno, dim, &
                                      & ele_loc(old_positions, 1), node_count(old_positions), &
                                      & ele_count(old_positions))

  end subroutine rtree_intersection_finder_set_input

  subroutine rtree_intersection_finder_find(new_positions, ele_B)
    type(vector_field), intent(in) :: new_positions
    integer, intent(in) :: ele_B

    integer :: dim, loc

    dim = new_positions%dim
    loc = ele_loc(new_positions, 1)

    call crtree_intersection_finder_find(reshape(ele_val(new_positions, ele_B), (/dim*loc/)), dim, loc)
    
  end subroutine rtree_intersection_finder_find
  
  function rtree_intersection_finder(positionsA, ndglnoA, &
       positionsB, ndglnoB, npredicates) result(map_ab)
    !!< As advancing_front_intersection_finder, but uses an rtree algorithm. For
    !!< testing *only*. For practical applications, use the linear algorithm.

    ! The positions and meshes of A and B
    real, dimension(:, :), intent(in) :: positionsA
    integer, dimension(:, :), intent(in) :: ndglnoA
    real, dimension(:, :), intent(in) :: positionsB
    integer, dimension(:, :), intent(in) :: ndglnoB
    integer, intent(out), optional :: npredicates
    ! for each element in A, the intersecting elements in B
    type(ilist), dimension(size(ndglnoA, 2)) :: map_ab

    integer :: i, j, id, nelms, ntests, dimA, verticesA, nodesA, dimB, nodesB, verticesB

    type(mesh_type) :: mesh
    type(vector_field), target :: lpositionsA, lpositionsB

    ewrite(1, *) "In rtree_intersection_finder"

    dimA = size(positionsA, 1)
    nodesA = size(positionsA, 2)
    verticesA = size(ndglnoA, 1)
    dimB = size(positionsB, 1)
    nodesB = size(positionsB, 2)
    verticesB = size(ndglnoB, 1)

    call allocate(mesh, dimA, nodesA, size(ndglnoA, 2), size(ndglnoA, 1))
    do i = 1, size(ndglnoA, 2)
      mesh%ndglno((i - 1) * verticesA + 1:i * verticesA) = ndglnoA(:, i)
    end do
    call allocate(lpositionsA, dimA, mesh)
    call deallocate(mesh)
    lpositionsA%val = positionsA
    
    call allocate(mesh, dimB, nodesB, size(ndglnoB, 2), size(ndglnoB, 1))
    do i = 1, size(ndglnoB, 2)
      mesh%ndglno((i - 1) * verticesB + 1:i * verticesB) = ndglnoB(:, i)
    end do
    call allocate(lpositionsB, dimB, mesh)
    call deallocate(mesh)
    lpositionsB%val = positionsB
    
!    call rtree_intersection_finder_set_input(positionsB, dimB, nodesB, size(ndglnoB, 1), verticesB, ndglnoB)
    call rtree_intersection_finder_set_input(lpositionsB)
    do i = 1, size(ndglnoA, 2)
      call rtree_intersection_finder_find(lpositionsA, i)
      call rtree_intersection_finder_query_output(nelms)
      do j = 1, nelms
        call rtree_intersection_finder_get_output(id, j)
        call insert(map_ab(i), id)
      end do
    end do
    call rtree_intersection_finder_reset(ntests)

    if (present(npredicates)) npredicates = ntests
    
    ewrite(1, *) "Exiting rtree_intersection_finder"
    
    call deallocate(lpositionsA)
    call deallocate(lpositionsB)
    
  end function rtree_intersection_finder
  
! IAKOVOS commented out
!  subroutine verify_map(mesh_field_a, mesh_field_b, map_ab, map_ab_reference)

  subroutine compute_bboxes(positionsB, bboxes_B)
    type(vector_field), intent(in) :: positionsB
    real, dimension(:, :, :), intent(out) :: bboxes_B
    integer :: ele_B

    do ele_B=1,ele_count(positionsB)
      bboxes_B(ele_B, :, :) = bbox(ele_val(positionsB, ele_B))
    end do
  end subroutine compute_bboxes

  function brute_force_search(posA, bboxes_B) result(map)
    real, dimension(:, :), intent(in) :: posA
    real, dimension(:, :, :), intent(in) :: bboxes_B

    type(ilist) :: map

    integer :: ele_B
    real, dimension(size(posA, 1), 2) :: bboxA

    bboxA = bbox(posA)

    do ele_B=1,size(bboxes_B, 1)
      if (bbox_predicate(bboxA, bboxes_B(ele_B, :, :))) then
        call insert(map, ele_B)
      end if
    end do

    if(map%length == 0) then
      FLAbort("Should never get here -- it has to intersect /something/!")
    end if
    
  end function brute_force_search

  function clueful_search(posA, possibles, bboxes_B, ele_A, neighbour) result(clues)
    real, dimension(:, :), intent(in) :: posA
    type(ilist), intent(in) :: possibles
    real, dimension(:, :, :), intent(in) :: bboxes_B
    integer, intent(in) :: ele_A, neighbour
    type(inode), pointer :: node
    type(ilist) :: clues
    real, dimension(size(posA, 1), 2) :: bboxA
    integer :: ele_B

    bboxA = bbox(posA)

    node => possibles%firstnode
    do while (associated(node))
      ele_B = node%value
      if (bbox_predicate(bboxA, bboxes_B(ele_B, :, :))) then
        call insert(clues, ele_B)
      end if
      node => node%next
    end do

  end function clueful_search

end module libsupermesh_intersection_finder_module
