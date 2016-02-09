#define BUF_SIZE 729
#include "fdebug.h"

module libsupermesh_tet_intersection_module

  use libsupermesh_debug
  use libsupermesh_tri_intersection_module

  implicit none

  private

  public :: tet_type, plane_type, intersect_tets_libwm, intersect_tets, &
    & intersect_tets_dt_public, intersect_tets_dt_real, get_planes

  type tet_type
    real, dimension(3, 4) :: V ! vertices of the tet
    integer, dimension(4) :: colours = -1 ! surface colours
  end type tet_type

  type plane_type
    real, dimension(3) :: normal
    real :: c
  end type plane_type

  type(tet_type), dimension(BUF_SIZE), save :: tet_array_tmp
  integer, save :: tet_cnt_tmp = 0
  logical, save :: mesh_allocated = .false.

  interface intersect_tets
    module procedure intersect_tets_dt, intersect_tets_dt_public, &
        intersect_tets_dt_real
  end interface intersect_tets

  interface get_planes
    module procedure get_planes_tet
  end interface get_planes

  integer, parameter, public :: tet_buf_size = BUF_SIZE

contains

  subroutine intersect_tets_libwm(tetA, tetB, nodesC, ndglnoC, n_tetsC)
    real, dimension(3, 4), intent(in) :: tetA
    real, dimension(3, 4), intent(in) :: tetB
    real, dimension(3, BUF_SIZE), intent(inout) :: nodesC
    integer, dimension(4, BUF_SIZE), intent(inout) :: ndglnoC
    integer, intent(out) :: n_tetsC

    integer :: i, nonods

    call cintersector_set_input(tetA, tetB, 3, 4)
    call cintersector_drive
    call cintersector_query(nonods, n_tetsC)
    call cintersector_get_output(nonods, n_tetsC, 3, 4, nodesC, ndglnoC)

  end subroutine intersect_tets_libwm

  subroutine intersect_tets_dt_real(tetA, tetB, tetsC, n_tetsC)
    real, dimension(3, 4), intent(in) :: tetA
    real, dimension(3, 4), intent(in) :: tetB
    real, dimension(3, 4, BUF_SIZE), intent(inout) :: tetsC
    integer, intent(out) :: n_tetsC

    integer :: i
    type(tet_type) :: tetA_t, tetB_t
    type(tet_type), dimension(BUF_SIZE), save :: tetsC_t

    tetA_t%v = tetA
    tetB_t%v = tetB
    call intersect_tets_dt_public(tetA_t, tetB_t, tetsC_t, n_tetsC)

    do i = 1, n_tetsC
      tetsC(:, :, i) = tetsC_t(i)%v
    end do

  end subroutine intersect_tets_dt_real

  subroutine intersect_tets_dt_public(tetA, tetB, tetsC, n_tetsC)
    type(tet_type), intent(in) :: tetA
    type(tet_type), intent(in) :: tetB
    type(tet_type), dimension(BUF_SIZE), intent(inout) :: tetsC
    integer, intent(out) :: n_tetsC

    call intersect_tets_dt(tetA, get_planes(tetB), tetsC, n_tetsC, volB = tet_volume(tetB))

  end subroutine intersect_tets_dt_public

  subroutine intersect_tets_dt(tetA, planesB, tetsC, n_tetsC, volB)
    type(tet_type), intent(in) :: tetA
    type(plane_type), dimension(4), intent(in)  :: planesB
    type(tet_type), dimension(BUF_SIZE), intent(inout) :: tetsC
    integer, intent(out) :: n_tetsC
    real, optional, intent(in) :: volB

    integer :: colour_tmp, i, j
    real, dimension(3) :: vec_tmp
    real :: tol, vol

    n_tetsC = 1
    tetsC(1) = tetA

    if(present(volB)) then
      tol = 10.0 * min(spacing(tet_volume(tetA)), spacing(volB))
    else
      tol = 10.0 * spacing(tet_volume(tetA))
    end if
    do i=1,size(planesB)
      ! Clip the tet_array against the i'th plane
      tet_cnt_tmp = 0

      do j=1,n_tetsC
        call clip(planesB(i), tetsC(j))
      end do

      if (i /= size(planesB)) then
        n_tetsC = tet_cnt_tmp
        tetsC(1:n_tetsC) = tet_array_tmp(1:n_tetsC)
      else
        ! Copy the result if the volume is >= tol
        n_tetsC = 0
        do j=1,tet_cnt_tmp
          vol = tet_volume(tet_array_tmp(j))

          if (vol >= tol) then
            n_tetsC = n_tetsC + 1
            tetsC(n_tetsC) = tet_array_tmp(j)
          end if
        end do
      end if
    end do

  end subroutine intersect_tets_dt

  subroutine clip(plane, tet)
  ! Clip tet against the plane
  ! and append any output to tet_array_tmp.
    type(plane_type), intent(in) :: plane
    type(tet_type), intent(in) :: tet

    real, dimension(4) :: dists
    integer :: neg_cnt, pos_cnt, zer_cnt
    integer, dimension(4) :: neg_idx, pos_idx, zer_idx
    integer :: i

    real :: invdiff, w0, w1
    type(tet_type) :: tet_tmp

    ! Negative == inside
    ! Positive == outside

    neg_cnt = 0
    pos_cnt = 0
    zer_cnt = 0

    dists = distances_to_plane(plane, tet)
    do i=1,4
      if (abs(dists(i)) < epsilon(0.0)) then
        zer_cnt = zer_cnt + 1
        zer_idx(zer_cnt) = i
      else if (dists(i) < 0.0) then
        neg_cnt = neg_cnt + 1
        neg_idx(neg_cnt) = i
      else if (dists(i) > 0.0) then
        pos_cnt = pos_cnt + 1
        pos_idx(pos_cnt) = i
      end if
    end do

    if (neg_cnt == 0) then
      ! tet is completely on positive side of plane, full clip
      return
    end if

    if (pos_cnt == 0) then
      ! tet is completely on negative side of plane, no clip
      tet_cnt_tmp = tet_cnt_tmp + 1
      tet_array_tmp(tet_cnt_tmp) = tet
      return
    end if

    ! The tet is split by the plane, so we have more work to do.

    select case(pos_cnt)
    case(3)
      ! +++-
      tet_cnt_tmp = tet_cnt_tmp + 1
      tet_array_tmp(tet_cnt_tmp) = tet
      do i=1,pos_cnt
        invdiff = 1.0 / ( dists(pos_idx(i)) - dists(neg_idx(1)) )
        w0 = -dists(neg_idx(1)) * invdiff
        w1 =  dists(pos_idx(i)) * invdiff
        tet_array_tmp(tet_cnt_tmp)%V(:, pos_idx(i)) = &
           w0 * tet_array_tmp(tet_cnt_tmp)%V(:, pos_idx(i)) + &
           w1 * tet_array_tmp(tet_cnt_tmp)%V(:, neg_idx(1))
      end do
      ! The colours will have been inherited already; we just need to zero
      ! the one corresponding to the plane cut
      tet_array_tmp(tet_cnt_tmp)%colours(face_no(pos_idx(1), pos_idx(2), pos_idx(3))) = 0
    case(2)
      select case(neg_cnt)
      case(2)
        ! ++--
        do i=1,pos_cnt
          invdiff = 1.0 / ( dists(pos_idx(i)) - dists(neg_idx(1)) )
          w0 = -dists(neg_idx(1)) * invdiff
          w1 =  dists(pos_idx(i)) * invdiff
          tet_tmp%V(:, i) = w0 * tet%V(:, pos_idx(i)) + w1 * tet%V(:, neg_idx(1))
        end do
        do i=1,neg_cnt
          invdiff = 1.0 / ( dists(pos_idx(i)) - dists(neg_idx(2)) )
          w0 = -dists(neg_idx(2)) * invdiff
          w1 =  dists(pos_idx(i)) * invdiff
          tet_tmp%V(:, i+2) = w0 * tet%V(:, pos_idx(i)) + w1 * tet%V(:, neg_idx(2))
        end do

        tet_cnt_tmp = tet_cnt_tmp + 1
        tet_array_tmp(tet_cnt_tmp) = tet
        tet_array_tmp(tet_cnt_tmp)%V(:, pos_idx(1)) = tet_tmp%V(:, 3)
        tet_array_tmp(tet_cnt_tmp)%V(:, pos_idx(2)) = tet_tmp%V(:, 2)
        tet_array_tmp(tet_cnt_tmp)%colours(neg_idx(1)) = 0
        tet_array_tmp(tet_cnt_tmp)%colours(neg_idx(2)) = 0

        tet_cnt_tmp = tet_cnt_tmp + 1
        tet_array_tmp(tet_cnt_tmp)%V(:, 1) = tet%V(:, neg_idx(2))
        tet_array_tmp(tet_cnt_tmp)%colours(1) = 0
        tet_array_tmp(tet_cnt_tmp)%V(:, 2) = tet_tmp%V(:, 4)
        tet_array_tmp(tet_cnt_tmp)%colours(2) = 0
        tet_array_tmp(tet_cnt_tmp)%V(:, 3) = tet_tmp%V(:, 3)
        tet_array_tmp(tet_cnt_tmp)%colours(3) = tet%colours(pos_idx(1))
        tet_array_tmp(tet_cnt_tmp)%V(:, 4) = tet_tmp%V(:, 2)
        tet_array_tmp(tet_cnt_tmp)%colours(4) = tet%colours(neg_idx(1))

        tet_cnt_tmp = tet_cnt_tmp + 1
        tet_array_tmp(tet_cnt_tmp)%V(:, 1) = tet%V(:, neg_idx(1))
        tet_array_tmp(tet_cnt_tmp)%colours(1) = 0
        tet_array_tmp(tet_cnt_tmp)%V(:, 2) = tet_tmp%V(:, 1)
        tet_array_tmp(tet_cnt_tmp)%colours(2) = 0
        tet_array_tmp(tet_cnt_tmp)%V(:, 3) = tet_tmp%V(:, 2)
        tet_array_tmp(tet_cnt_tmp)%colours(3) = tet%colours(pos_idx(2))
        tet_array_tmp(tet_cnt_tmp)%V(:, 4) = tet_tmp%V(:, 3)
        tet_array_tmp(tet_cnt_tmp)%colours(4) = tet%colours(neg_idx(2))
      case(1)
        ! ++-0
        tet_cnt_tmp = tet_cnt_tmp + 1
        tet_array_tmp(tet_cnt_tmp) = tet
        do i=1,pos_cnt
          invdiff = 1.0 / ( dists(pos_idx(i)) - dists(neg_idx(1)) )
          w0 = -dists(neg_idx(1)) * invdiff
          w1 =  dists(pos_idx(i)) * invdiff
          tet_array_tmp(tet_cnt_tmp)%V(:, pos_idx(i)) = &
             w0 * tet_array_tmp(tet_cnt_tmp)%V(:, pos_idx(i)) + &
             w1 * tet_array_tmp(tet_cnt_tmp)%V(:, neg_idx(1))
        end do
        tet_array_tmp(tet_cnt_tmp)%colours(neg_idx(1)) = 0
      end select
    case(1)
      select case(neg_cnt)
      case(3)
        ! +---
        do i=1,neg_cnt
          invdiff = 1.0 / ( dists(pos_idx(1)) - dists(neg_idx(i)) )
          w0 = -dists(neg_idx(i)) * invdiff
          w1 =  dists(pos_idx(1)) * invdiff
          tet_tmp%V(:, i) = w0 * tet%V(:, pos_idx(1)) + w1 * tet%V(:, neg_idx(i))
        end do

        tet_cnt_tmp = tet_cnt_tmp + 1
        tet_array_tmp(tet_cnt_tmp) = tet
        tet_array_tmp(tet_cnt_tmp)%V(:, pos_idx(1)) = tet_tmp%V(:, 1)
        tet_array_tmp(tet_cnt_tmp)%colours(neg_idx(1)) = 0

        tet_cnt_tmp = tet_cnt_tmp + 1
        tet_array_tmp(tet_cnt_tmp)%V(:, 1) = tet_tmp%V(:, 1)
        tet_array_tmp(tet_cnt_tmp)%colours(1) = tet%colours(neg_idx(1))
        tet_array_tmp(tet_cnt_tmp)%V(:, 2) = tet%V(:, neg_idx(2))
        tet_array_tmp(tet_cnt_tmp)%colours(2) = 0
        tet_array_tmp(tet_cnt_tmp)%V(:, 3) = tet%V(:, neg_idx(3))
        tet_array_tmp(tet_cnt_tmp)%colours(3) = tet%colours(neg_idx(3))
        tet_array_tmp(tet_cnt_tmp)%V(:, 4) = tet_tmp%V(:, 2)
        tet_array_tmp(tet_cnt_tmp)%colours(4) = 0

        tet_cnt_tmp = tet_cnt_tmp + 1
        tet_array_tmp(tet_cnt_tmp)%V(:, 1) = tet%V(:, neg_idx(3))
        tet_array_tmp(tet_cnt_tmp)%colours(1) = 0
        tet_array_tmp(tet_cnt_tmp)%V(:, 2) = tet_tmp%V(:, 2)
        tet_array_tmp(tet_cnt_tmp)%colours(2) = tet%colours(neg_idx(2))
        tet_array_tmp(tet_cnt_tmp)%V(:, 3) = tet_tmp%V(:, 3)
        tet_array_tmp(tet_cnt_tmp)%colours(3) = 0
        tet_array_tmp(tet_cnt_tmp)%V(:, 4) = tet_tmp%V(:, 1)
        tet_array_tmp(tet_cnt_tmp)%colours(4) = tet%colours(neg_idx(1))
      case(2)
        ! +--0
        do i=1,neg_cnt
          invdiff = 1.0 / ( dists(pos_idx(1)) - dists(neg_idx(i)) )
          w0 = -dists(neg_idx(i)) * invdiff
          w1 =  dists(pos_idx(1)) * invdiff
          tet_tmp%V(:, i) = w0 * tet%V(:, pos_idx(1)) + w1 * tet%V(:, neg_idx(i))
        end do

        tet_cnt_tmp = tet_cnt_tmp + 1
        tet_array_tmp(tet_cnt_tmp) = tet
        tet_array_tmp(tet_cnt_tmp)%V(:, pos_idx(1)) = tet_tmp%V(:, 1)
        tet_array_tmp(tet_cnt_tmp)%colours(neg_idx(1)) = 0

        tet_cnt_tmp = tet_cnt_tmp + 1
        tet_array_tmp(tet_cnt_tmp)%V(:, 1) = tet_tmp%V(:, 2)
        tet_array_tmp(tet_cnt_tmp)%colours(1) = 0
        tet_array_tmp(tet_cnt_tmp)%V(:, 2) = tet%V(:, zer_idx(1))
        tet_array_tmp(tet_cnt_tmp)%colours(2) = tet%colours(zer_idx(1))
        tet_array_tmp(tet_cnt_tmp)%V(:, 3) = tet%V(:, neg_idx(2))
        tet_array_tmp(tet_cnt_tmp)%colours(3) = 0
        tet_array_tmp(tet_cnt_tmp)%V(:, 4) = tet_tmp%V(:, 1)
        tet_array_tmp(tet_cnt_tmp)%colours(4) = tet%colours(neg_idx(1))
      case(1)
        ! +-00
        invdiff = 1.0 / ( dists(pos_idx(1)) - dists(neg_idx(1)) )
        w0 = -dists(neg_idx(1)) * invdiff
        w1 =  dists(pos_idx(1)) * invdiff

        tet_cnt_tmp = tet_cnt_tmp + 1
        tet_array_tmp(tet_cnt_tmp) = tet
        tet_array_tmp(tet_cnt_tmp)%V(:, pos_idx(1)) = w0 * tet%V(:, pos_idx(1)) + w1 * tet%V(:, neg_idx(1))
        tet_array_tmp(tet_cnt_tmp)%colours(neg_idx(1)) = 0
      end select
    end select

  end subroutine clip
  
  pure function get_planes_tet(tet) result(plane)
    type(tet_type), intent(in) :: tet
    type(plane_type), dimension(4) :: plane

    real, dimension(3) :: edge10, edge20, edge30, edge21, edge31
    real :: det
    integer :: i

    edge10 = tet%V(:, 2) - tet%V(:, 1);
    edge20 = tet%V(:, 3) - tet%V(:, 1);
    edge30 = tet%V(:, 4) - tet%V(:, 1);
    edge21 = tet%V(:, 3) - tet%V(:, 2);
    edge31 = tet%V(:, 4) - tet%V(:, 2);

    plane(1)%normal = unit_cross(edge20, edge10)
    plane(2)%normal = unit_cross(edge10, edge30)
    plane(3)%normal = unit_cross(edge30, edge20)
    plane(4)%normal = unit_cross(edge21, edge31)

    det = dot_product(edge10, plane(4)%normal)
    if (det < 0) then
      do i=1,4
        plane(i)%normal = -plane(i)%normal
      end do
    end if

    ! And calibrate what is the zero of this plane by dotting with
    ! a point we know to be on it
    do i=1,4
      plane(i)%c = dot_product(tet%V(:, i), plane(i)%normal)
    end do

  end function get_planes_tet

!   function get_planes_hex(positions, ele) result(plane)
!     type(vector_field), intent(in) :: positions
!     integer, intent(in) :: ele
!     type(plane_type), dimension(6) :: plane
!     integer, dimension(:), pointer :: faces
!     integer :: i, face
!     integer, dimension(4) :: fnodes
!     real, dimension(positions%dim, face_ngi(positions, ele)) :: normals
! 
!     ! This could be done much more efficiently by exploiting
!     ! more information about how we number faces and such on a hex
! 
!     assert(positions%mesh%shape%numbering%family == FAMILY_CUBE)
!     assert(positions%mesh%faces%shape%numbering%family == FAMILY_CUBE)
!     assert(positions%mesh%shape%degree == 1)
!     assert(has_faces(positions%mesh))
! 
!     faces => ele_faces(positions, ele)
!     assert(size(faces) == 6)
! 
!     do i=1,size(faces)
!       face = faces(i)
!       fnodes = face_global_nodes(positions, face)
! 
!       call transform_facet_to_physical(positions, face, normal=normals)
!       plane(i)%normal = normals(:, 1)
! 
!       ! Now we calibrate the constant (setting the 'zero level' of the plane, as it were)
!       ! with a node we know is on the face
!       plane(i)%c = dot_product(plane(i)%normal, node_val(positions, fnodes(1)))
! 
!     end do
!   end function get_planes_hex

  pure function unit_cross(vecA, vecB) result(cross)
    real, dimension(3), intent(in) :: vecA, vecB
    real, dimension(3) :: cross
    cross(1) = vecA(2) * vecB(3) - vecA(3) * vecB(2)
    cross(2) = vecA(3) * vecB(1) - vecA(1) * vecB(3)
    cross(3) = vecA(1) * vecB(2) - vecA(2) * vecB(1)

    cross = cross / norm2(cross)
  end function unit_cross

  pure function distances_to_plane(plane, tet) result(dists)
    type(plane_type), intent(in) :: plane
    type(tet_type), intent(in) :: tet
    real, dimension(4) :: dists
    integer :: i

    forall(i=1:4)
      dists(i) = dot_product(plane%normal, tet%V(:, i)) - plane%c
    end forall
  end function distances_to_plane

  pure function tet_volume(tet) result(vol)
    type(tet_type), intent(in) :: tet
    real :: vol
    real, dimension(3) :: cross, vecA, vecB, vecC

    vecA = tet%V(:, 1) - tet%V(:, 4)
    vecB = tet%V(:, 2) - tet%V(:, 4)
    vecC = tet%V(:, 3) - tet%V(:, 4)

    cross(1) = vecB(2) * vecC(3) - vecB(3) * vecC(2)
    cross(2) = vecB(3) * vecC(1) - vecB(1) * vecC(3)
    cross(3) = vecB(1) * vecC(2) - vecB(2) * vecC(1)

    vol = abs(dot_product(vecA, cross)) / 6.0
  end function tet_volume

  function face_no(i, j, k) result(face)
    ! Given three local node numbers, what is the face that they share?
    integer, intent(in) :: i, j, k
    integer :: face

    do face=1,4
      if (face /= i .and. face /= j .and. face /= k) return
    end do

  end function face_no

end module libsupermesh_tet_intersection_module
