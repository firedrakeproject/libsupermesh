!#define COUNT_INTERSECTION_TESTS
#include "fdebug.h"

module libsupermesh_intersection_finder

  use iso_c_binding
  use libsupermesh_fields_dummy
  use libsupermesh_linked_lists

  implicit none

  private

  interface crtree_intersection_finder_set_input
    subroutine libsupermesh_cintersection_finder_set_input(positions, enlist, dim, loc, nnodes, nelements) bind(c)
      use iso_c_binding
      implicit none
      integer(kind = c_int), intent(in) :: dim, loc, nnodes, nelements
      real(kind = c_double), intent(in), dimension(dim, nnodes) :: positions
      integer(kind = c_int), intent(in), dimension(loc, nelements) :: enlist
    end subroutine libsupermesh_cintersection_finder_set_input
  end interface crtree_intersection_finder_set_input

  interface crtree_intersection_finder_find
    subroutine libsupermesh_cintersection_finder_find(positions, dim, loc) bind(c)
      use iso_c_binding
      implicit none
      integer(kind = c_int), intent(in) :: dim, loc
      real(kind = c_double), dimension(dim, loc) :: positions
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

  public :: rtree_intersection_finder_set_input, &
    & rtree_intersection_finder_find, rtree_intersection_finder_query_output, &
    & rtree_intersection_finder_get_output, rtree_intersection_finder_reset

!   interface
!     function mpi_wtime()
!       implicit none
!       real :: mpi_wtime
!     end function mpi_wtime
!   end interface
  
#ifdef COUNT_INTERSECTION_TESTS
  integer, save :: intersection_tests_count = 0
#endif
  public :: intersection_tests, reset_intersection_tests_counter

  type intersections
    integer, dimension(:), pointer :: v
    integer :: n
  end type intersections

  interface deallocate
    module procedure deallocate_intersections, deallocate_intersections_vector
  end interface deallocate
  
  public :: intersections, ilist, inode, deallocate

  interface intersection_finder
    module procedure advancing_front_intersection_finder_intersections, &
      & advancing_front_intersection_finder_csr_sparsity, &
      & advancing_front_intersection_finder_lists
  end interface intersection_finder

  interface advancing_front_intersection_finder
    module procedure advancing_front_intersection_finder_intersections, &
      & advancing_front_intersection_finder_csr_sparsity, &
      & advancing_front_intersection_finder_lists
  end interface advancing_front_intersection_finder

  interface rtree_intersection_finder
    module procedure rtree_intersection_finder_intersections, &
      & rtree_intersection_finder_csr_sparsity, &
      & rtree_intersection_finder_lists
  end interface rtree_intersection_finder

  interface brute_force_intersection_finder
    module procedure brute_force_intersection_finder_intersections, &
      & brute_force_intersection_finder_csr_sparsity, &
      & brute_force_intersection_finder_lists
  end interface brute_force_intersection_finder

  public :: intersection_finder, advancing_front_intersection_finder, &
    & rtree_intersection_finder, brute_force_intersection_finder

contains

  function intersection_tests() result(tests)
    !!< Return the number of intersection tests

    integer :: tests

#ifdef COUNT_INTERSECTION_TESTS
    tests = intersection_tests_count
#else
    FLAbort("Counting of intersection tests is not available")
    ! To keep the compiler quiet
    tests = -1
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

  pure subroutine deallocate_intersections(ints)
    type(intersections), intent(inout) :: ints

    deallocate(ints%v)

  end subroutine deallocate_intersections

  pure subroutine deallocate_intersections_vector(ints)
    type(intersections), dimension(:), intent(inout) :: ints

    integer :: i

    do i = 1, size(ints)
      deallocate(ints(i)%v)
    end do

  end subroutine deallocate_intersections_vector

  pure function bbox(coords)
    ! dim x loc
    real, dimension(:, :), intent(in) :: coords

    real, dimension(2, size(coords, 1)) :: bbox

    integer :: i, j

    bbox(1, :) = coords(:, 1)
    bbox(2, :) = coords(:, 1)
    do i = 2, size(coords, 2)
      do j = 1, size(coords, 1)
        bbox(1, j) = min(bbox(1, j), coords(j, i))
        bbox(2, j) = max(bbox(2, j), coords(j, i))
      end do
    end do

  end function bbox

  pure function bboxes_intersect(bbox_1, bbox_2) result(intersect)
    ! 2 x dim
    real, dimension(:, :), intent(in) :: bbox_1
    ! 2 x dim
    real, dimension(:, :), intent(in) :: bbox_2

    logical :: intersect

    integer :: i

#ifdef COUNT_INTERSECTION_TESTS
    intersection_tests_count = intersection_tests_count + 1
#endif

    do i = 1, size(bbox_1, 2)
      ! Strict inequalities required here for the advancing front intersection
      ! finder to work with axis aligned elements
      if(bbox_2(2, i) < bbox_1(1, i) .or. bbox_2(1, i) > bbox_1(2, i)) then
        intersect = .false.
        return
      end if
    end do
    intersect = .true.

  end function bboxes_intersect

  ! Advancing front intersection finder, as described in P. E. Farrell and
  ! J. R. Maddison, "Conservative interpolation between volume meshes by local
  ! Galerkin projection", Computer Methods in Applied Mechanics and Engineering,
  ! 200, pp. 89--100, 2011
  subroutine advancing_front_intersection_finder_intersections(positions_a, enlist_a, positions_b, enlist_b, map_ab)
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

    integer :: i, j
    integer :: dim, nnodes_b, nnodes_a, nelements_b, nelements_a

    real, dimension(2, size(positions_a, 1)) :: bbox_a
    real, dimension(:, :, :), allocatable :: bboxes_b
    type(eelist_type) :: eelist_b, eelist_a

    integer :: clue_a, ele_b, ele_a, loc_b, loc_a, neigh_b, neigh_a, nsub_a, seed_a
    integer, dimension(:), allocatable :: front_b, ints
    logical, dimension(:), allocatable :: seen_b, seen_a
    integer, dimension(:, :), allocatable :: front_a
    integer :: nfront_b, nfront_a, nints

!     real :: t_0

    dim = size(positions_b, 1)
    nnodes_b = size(positions_b, 2)
    nnodes_a = size(positions_a, 2)
    loc_b = size(enlist_b, 1)
    nelements_b = size(enlist_b, 2)
    loc_a = size(enlist_a, 1)
    nelements_a = size(enlist_a, 2)

!     t_0 = mpi_wtime()
    call mesh_eelist(nnodes_b, enlist_b, sloc(dim, loc_b), eelist_b)
    call mesh_eelist(nnodes_a, enlist_a, sloc(dim, loc_a), eelist_a)
!     ewrite(2, "(a,e25.17e3)") "eelist creation time = ", mpi_wtime() - t_0
    allocate(bboxes_b(2, dim, nelements_b))
    do i = 1, nelements_b
      bboxes_b(:, :, i) = bbox(positions_b(:, enlist_b(:, i)))
    end do

    allocate(seen_b(nelements_b), seen_a(nelements_a), front_b(nelements_b), &
      & front_a(2, nelements_a), ints(nelements_b))
    seen_b = .false.
    seen_a = .false.
    ! Stage 0: Initial target mesh seed element
    seed_a = 1
    nsub_a = 0
    seed_a_loop: do
      nsub_a = nsub_a + 1
      seen_a(seed_a) = .true.
      bbox_a = bbox(positions_a(:, enlist_a(:, seed_a)))

      ! Stage 1a: Find intersections with the target mesh seed element via a
      ! brute force search
      nints = 0
      do ele_b = 1, nelements_b
        if(bboxes_intersect(bbox_a, bboxes_b(:, :, ele_b))) then
          nints = nints + 1
          ints(nints) = ele_b
        end if
      end do
      allocate(map_ab(seed_a)%v(nints))
      map_ab(seed_a)%v = ints(:nints)
      map_ab(seed_a)%n = nints

      ! Stage 1b: Advance the target mesh front
      nfront_a = 0
      do i = 1, eelist_a%n(seed_a)
        neigh_a = eelist_a%v(i, seed_a)
        if(.not. seen_a(neigh_a)) then
          nfront_a = nfront_a + 1
          front_a(1, nfront_a) = neigh_a
          front_a(2, nfront_a) = seed_a
          seen_a(neigh_a) = .true.
        end if
      end do

      do while(nfront_a > 0)
        ele_a = front_a(1, nfront_a)
        clue_a = front_a(2, nfront_a)
        nfront_a = nfront_a - 1
        bbox_a = bbox(positions_a(:, enlist_a(:, ele_a)))

        ! Stage 2a: Initialise the donor mesh front
        nfront_b = map_ab(clue_a)%n
        front_b(:nfront_b) = map_ab(clue_a)%v
        seen_b(front_b(:nfront_b)) = .true.

        ! Stage 2b: Find intersections with the target mesh element by
        ! advancing the donor mesh front
        nints = 0
        i = 1
        do while(i <= nfront_b)
          ele_b = front_b(i)
          if(bboxes_intersect(bbox_a, bboxes_b(:, :, ele_b))) then
            ! An intersection has been found
            nints = nints + 1
            ints(nints) = ele_b
            ! Advance the donor mesh front
            do j = 1, eelist_b%n(ele_b)
              neigh_b = eelist_b%v(j, ele_b)
              if(.not. seen_b(neigh_b)) then
                nfront_b = nfront_b + 1
                front_b(nfront_b) = neigh_b
                seen_b(neigh_b) = .true.
              end if
            end do
          end if
          i = i + 1
        end do
        do i = 1, nfront_b
          seen_b(front_b(i)) = .false.
        end do
!         if(nints == 0) then
!           ewrite(-1, "(a,i0)") "WARNING: Failed to find intersections for target element ", ele_a
!         end if
        allocate(map_ab(ele_a)%v(nints))
        map_ab(ele_a)%v = ints(:nints)
        map_ab(ele_a)%n = nints

        ! Stage 2c: Advance the target mesh front
        do i = 1, eelist_a%n(ele_a)
          neigh_a = eelist_a%v(i, ele_a)
          if(.not. seen_a(neigh_a)) then
            nfront_a = nfront_a + 1
            front_a(1, nfront_a) = neigh_a
            front_a(2, nfront_a) = ele_a
            seen_a(neigh_a) = .true.
          end if
        end do
      end do

      ! Stage 3: Find a new target mesh seed
      do while(seen_a(seed_a))
        seed_a = seed_a + 1
        if(seed_a > nelements_a) exit seed_a_loop
      end do
      if(nsub_a == 1) then
        ewrite(-1, "(a)") "WARNING: Target mesh is not connected"
      end if
    end do seed_a_loop
    if(nsub_a > 1) then
      ewrite(-1 ,"(a,i0)") "WARNING: Number of target connected sub-domains = ", nsub_a
    end if

    deallocate(seen_b, seen_a, front_b, front_a, ints)
    call deallocate(eelist_b)
    call deallocate(eelist_a)
    deallocate(bboxes_b)

  end subroutine advancing_front_intersection_finder_intersections

  pure subroutine intersections_to_csr_sparsity(ints, indices, indptr)
    type(intersections), dimension(:), intent(in) :: ints
    integer, dimension(:), allocatable, intent(out) :: indices
    integer, dimension(:), intent(out) :: indptr

    integer :: i, j, n

    n = 0
    do i = 1, size(ints)
      n = n + ints(i)%n
    end do

    allocate(indices(n))
    indptr(1) = 1
    do i = 1, size(ints)
      indptr(i + 1) = indptr(i) + ints(i)%n
      do j = 1, ints(i)%n
        indices(indptr(i) + j - 1) = ints(i)%v(j)
      end do
    end do

  end subroutine intersections_to_csr_sparsity

  subroutine advancing_front_intersection_finder_csr_sparsity(positions_a, enlist_a, positions_b, enlist_b, map_ab_indices, map_ab_indptr)
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

    type(intersections), dimension(:), allocatable :: map_ab

    allocate(map_ab(size(enlist_a, 2)))
    call advancing_front_intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
    call intersections_to_csr_sparsity(map_ab, map_ab_indices, map_ab_indptr)
    call deallocate(map_ab)
    deallocate(map_ab)

  end subroutine advancing_front_intersection_finder_csr_sparsity

  subroutine intersections_to_lists(ints, lists)
    type(intersections), dimension(:), intent(in) :: ints
    type(ilist), dimension(:), intent(out) :: lists

    integer :: i, j

    do i = 1, size(ints)
      do j = 1, ints(i)%n
        call insert(lists(i), ints(i)%v(j))
      end do
    end do

  end subroutine intersections_to_lists

  subroutine advancing_front_intersection_finder_lists(positions_a, enlist_a, positions_b, enlist_b, map_ab)
    ! dim x nnodes_a
    real, dimension(:, :), intent(in) :: positions_a
    ! loc_a x nelements_a
    integer, dimension(:, :), intent(in) :: enlist_a
    ! dim x nnodes_b
    real, dimension(:, :), intent(in) :: positions_b
    ! loc_b x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    ! nelements_a
    type(ilist), dimension(:), intent(out) :: map_ab

    type(intersections), dimension(:), allocatable :: lmap_ab

    allocate(lmap_ab(size(enlist_a, 2)))
    call advancing_front_intersection_finder(positions_a, enlist_a, positions_b, enlist_b, lmap_ab)
    call intersections_to_lists(lmap_ab, map_ab)
    call deallocate(lmap_ab)
    deallocate(lmap_ab)

  end subroutine advancing_front_intersection_finder_lists

  subroutine rtree_intersection_finder_set_input(positions, enlist)
    ! dim x nnodes_b
    real(kind = c_double), dimension(:, :), intent(in) :: positions
    ! loc x nelements_b
    integer(kind = c_int), dimension(:, :), intent(in) :: enlist

    integer(kind = c_int) :: loc, dim, nelements, nnodes

    dim = size(positions, 1)
    nnodes = size(positions, 2)
    loc = size(enlist, 1)
    nelements = size(enlist, 2)

    call crtree_intersection_finder_set_input(positions, enlist, dim, loc, nnodes, nelements)

  end subroutine rtree_intersection_finder_set_input

  subroutine rtree_intersection_finder_find(positions)
    ! dim x loc
    real(kind = c_double), dimension(:, :), intent(in) :: positions

    integer(kind = c_int) :: loc, dim

    dim = size(positions, 1)
    loc = size(positions, 2)

    call crtree_intersection_finder_find(positions, dim, loc)

  end subroutine rtree_intersection_finder_find

  subroutine rtree_intersection_finder_intersections(positions_a, enlist_a, positions_b, enlist_b, map_ab)
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

    integer :: ele_a, ele_b, i, nelements_a, nints, ntests

    nelements_a = size(enlist_a, 2)

    call rtree_intersection_finder_set_input(positions_b, enlist_b)
    do ele_a = 1, nelements_a
      call rtree_intersection_finder_find(positions_a(:, enlist_a(:, ele_a)))
      call rtree_intersection_finder_query_output(nints)
      map_ab(ele_a)%n = nints
      allocate(map_ab(ele_a)%v(nints))
      do i = 1, nints
        call rtree_intersection_finder_get_output(ele_b, i)
        map_ab(ele_a)%v(i) = ele_b
      end do
    end do
    call rtree_intersection_finder_reset(ntests)
#ifdef COUNT_INTERSECTION_TESTS
    intersection_tests_count = intersection_tests_count + ntests
#endif

  end subroutine rtree_intersection_finder_intersections

  subroutine rtree_intersection_finder_csr_sparsity(positions_a, enlist_a, positions_b, enlist_b, map_ab_indices, map_ab_indptr)
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

    type(intersections), dimension(:), allocatable :: map_ab

    allocate(map_ab(size(enlist_a, 2)))
    call rtree_intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
    call intersections_to_csr_sparsity(map_ab, map_ab_indices, map_ab_indptr)
    call deallocate(map_ab)
    deallocate(map_ab)

  end subroutine rtree_intersection_finder_csr_sparsity

  subroutine rtree_intersection_finder_lists(positions_a, enlist_a, positions_b, enlist_b, map_ab)
    ! dim x nnodes_a
    real, dimension(:, :), intent(in) :: positions_a
    ! loc_a x nelements_a
    integer, dimension(:, :), intent(in) :: enlist_a
    ! dim x nnodes_b
    real, dimension(:, :), intent(in) :: positions_b
    ! loc_b x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    ! nelements_a
    type(ilist), dimension(:), intent(out) :: map_ab

    type(intersections), dimension(:), allocatable :: lmap_ab

    allocate(lmap_ab(size(enlist_a, 2)))
    call rtree_intersection_finder(positions_a, enlist_a, positions_b, enlist_b, lmap_ab)
    call intersections_to_lists(lmap_ab, map_ab)
    call deallocate(lmap_ab)
    deallocate(lmap_ab)

  end subroutine rtree_intersection_finder_lists

  ! Brute force intersection finder.
  pure subroutine brute_force_intersection_finder_intersections(positions_a, enlist_a, positions_b, enlist_b, map_ab)
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

    integer :: ele_a, ele_b, nelements_a, nelements_b
    real, dimension(2, size(positions_a, 1)) :: bbox_a, bbox_b

    integer, dimension(:), allocatable :: ints
    integer :: nints

    nelements_a = size(enlist_a, 2)
    nelements_b = size(enlist_b, 2)

    allocate(ints(nelements_b))
    do ele_a = 1, nelements_a
      bbox_a = bbox(positions_a(:, enlist_a(:, ele_a)))
      nints = 0
      do ele_b = 1, nelements_b
        bbox_b = bbox(positions_b(:, enlist_b(:, ele_b)))
        if(bboxes_intersect(bbox_a, bbox_b)) then
          nints = nints + 1
          ints(nints) = ele_b
        end if
      end do

      map_ab(ele_a)%n = nints
      allocate(map_ab(ele_a)%v(nints))
      map_ab(ele_a)%v = ints(:nints)
    end do
    deallocate(ints)

  end subroutine brute_force_intersection_finder_intersections

  pure subroutine brute_force_intersection_finder_csr_sparsity(positions_a, enlist_a, positions_b, enlist_b, map_ab_indices, map_ab_indptr)
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

    type(intersections), dimension(:), allocatable :: map_ab

    allocate(map_ab(size(enlist_a, 2)))
    call brute_force_intersection_finder(positions_a, enlist_a, positions_b, enlist_b, map_ab)
    call intersections_to_csr_sparsity(map_ab, map_ab_indices, map_ab_indptr)
    call deallocate(map_ab)
    deallocate(map_ab)

  end subroutine brute_force_intersection_finder_csr_sparsity

  subroutine brute_force_intersection_finder_lists(positions_a, enlist_a, positions_b, enlist_b, map_ab)
    ! dim x nnodes_a
    real, dimension(:, :), intent(in) :: positions_a
    ! loc_a x nelements_a
    integer, dimension(:, :), intent(in) :: enlist_a
    ! dim x nnodes_b
    real, dimension(:, :), intent(in) :: positions_b
    ! loc_b x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    ! nelements_a
    type(ilist), dimension(:), intent(out) :: map_ab

    type(intersections), dimension(:), allocatable :: lmap_ab

    allocate(lmap_ab(size(enlist_a, 2)))
    call brute_force_intersection_finder(positions_a, enlist_a, positions_b, enlist_b, lmap_ab)
    call intersections_to_lists(lmap_ab, map_ab)
    call deallocate(lmap_ab)
    deallocate(lmap_ab)

  end subroutine brute_force_intersection_finder_lists

end module libsupermesh_intersection_finder