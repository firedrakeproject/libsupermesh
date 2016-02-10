#include "fdebug.h"

#define BUF_SIZE 81

module libsupermesh_tet_intersection

  implicit none

  private

  public :: tet_type, plane_type, max_n_tetsC, intersect_tets, get_planes, &
    & tetrahedron_volume

  type tet_type
    real, dimension(3, 4) :: v ! vertices of the tet
    integer, dimension(4) :: colours = -1 ! surface colours
  end type tet_type

  type plane_type
    real, dimension(3) :: normal
    real :: c
  end type plane_type

  interface intersect_tets
    module procedure intersect_tets_real, intersect_tets_tet, &
      & intersect_tets_planes
  end interface intersect_tets

  interface get_planes
    module procedure get_planes_tet
  end interface get_planes
  
  interface tetrahedron_volume
    module procedure tetrahedron_volume_real, tetrahedron_volume_tet
  end interface tetrahedron_volume

  integer, parameter, public :: tet_buf_size = BUF_SIZE

  type(tet_type), dimension(BUF_SIZE), save :: tets_tmp_buf
  integer, save :: n_tets_tmp = 0

contains

  subroutine intersect_tets_real(tetA, tetB, tetsC, n_tetsC)
    real, dimension(3, 4), intent(in) :: tetA
    real, dimension(3, 4), intent(in) :: tetB
    real, dimension(3, 4, BUF_SIZE), intent(inout) :: tetsC
    integer, intent(out) :: n_tetsC

    integer :: i
    type(tet_type) :: tetA_t, tetB_t
    type(tet_type), dimension(BUF_SIZE), save :: tetsC_t

    tetA_t%v = tetA
    tetB_t%v = tetB
    call intersect_tets(tetA_t, tetB_t, tetsC_t, n_tetsC)

    do i = 1, n_tetsC
      tetsC(:, :, i) = tetsC_t(i)%v
    end do

  end subroutine intersect_tets_real

  subroutine intersect_tets_tet(tetA, tetB, tetsC, n_tetsC)
    type(tet_type), intent(in) :: tetA
    type(tet_type), intent(in) :: tetB
    type(tet_type), dimension(BUF_SIZE), intent(inout) :: tetsC
    integer, intent(out) :: n_tetsC

    call intersect_tets_planes_buf(tetA, get_planes(tetB), tetsC, n_tetsC, volB = tetrahedron_volume(tetB))

  end subroutine intersect_tets_tet

  subroutine intersect_tets_planes_buf(tetA, planesB, tetsC, n_tetsC, volB)
    type(tet_type), intent(in) :: tetA
    type(plane_type), dimension(4), intent(in)  :: planesB
    type(tet_type), dimension(BUF_SIZE), intent(inout) :: tetsC
    integer, intent(out) :: n_tetsC
    real, optional, intent(in) :: volB

    integer :: i, j
    real :: tol, vol

    n_tetsC = 1
    tetsC(1) = tetA

    if(present(volB)) then
      tol = 10.0 * min(spacing(tetrahedron_volume(tetA)), spacing(volB))
    else
      tol = 10.0 * spacing(tetrahedron_volume(tetA))
    end if
    do i = 1, size(planesB)
      ! Clip the tet_array against the i'th plane
      n_tets_tmp = 0

      do j = 1, n_tetsC
        call clip_buf(planesB(i), tetsC(j))
      end do

      if(i /= size(planesB)) then
        n_tetsC = n_tets_tmp
        tetsC(:n_tetsC) = tets_tmp_buf(:n_tetsC)
      else
        ! Copy the result if the volume is >= tol
        n_tetsC = 0
        do j = 1, n_tets_tmp
          vol = tetrahedron_volume(tets_tmp_buf(j))
          if(vol >= tol) then
            n_tetsC = n_tetsC + 1
            tetsC(n_tetsC) = tets_tmp_buf(j)
          end if
        end do
      end if
    end do

  end subroutine intersect_tets_planes_buf
  
  pure function max_n_tetsC(n_planesB)
    integer, intent(in) :: n_planesB
    
    integer :: max_n_tetsC
    
    max_n_tetsC = 3 ** n_planesB
    
  end function max_n_tetsC
  
  subroutine intersect_tets_planes(tetA, planesB, tetsC, n_tetsC, volB, work)
    type(tet_type), intent(in) :: tetA
    type(plane_type), dimension(:), intent(in)  :: planesB
    type(tet_type), dimension(:), intent(inout) :: tetsC
    integer, intent(out) :: n_tetsC
    real, optional, intent(in) :: volB
    type(tet_type), dimension(:), target, optional, intent(inout) :: work

    integer :: i, j, ntets_new
    real :: tol, vol
    type(tet_type), dimension(:), pointer :: tets_new
    
    if(present(work)) then
      tets_new => work
    else
      allocate(tets_new(max_n_tetsC(size(planesB))))
    end if

    n_tetsC = 1
    tetsC(1) = tetA

    if(present(volB)) then
      tol = 10.0 * min(spacing(tetrahedron_volume(tetA)), spacing(volB))
    else
      tol = 10.0 * spacing(tetrahedron_volume(tetA))
    end if
    do i = 1, size(planesB)
      ! Clip the tet_array against the i'th plane
      ntets_new = 0

      do j = 1, n_tetsC
        call clip(planesB(i), tetsC(j), tets_new, ntets_new)
      end do

      if(i /= size(planesB)) then
        n_tetsC = ntets_new
        tetsC(:n_tetsC) = tets_new(:n_tetsC)
      else
        ! Copy the result if the volume is >= tol
        n_tetsC = 0
        do j = 1, ntets_new
          vol = tetrahedron_volume(tets_new(j))
          if(vol >= tol) then
            n_tetsC = n_tetsC + 1
            tetsC(n_tetsC) = tets_new(j)
          end if
        end do
      end if
    end do
    
    if(.not. present(work)) deallocate(tets_new)
  
  end subroutine intersect_tets_planes

  subroutine clip_buf(plane, tet)
    ! Clip tet against the plane and append any output to tets_tmp_buf.
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
      n_tets_tmp = n_tets_tmp + 1
      tets_tmp_buf(n_tets_tmp) = tet
      return
    end if

    ! The tet is split by the plane, so we have more work to do.

    select case(pos_cnt)
    case(3)
      ! +++-
      n_tets_tmp = n_tets_tmp + 1
      tets_tmp_buf(n_tets_tmp) = tet
      do i=1,pos_cnt
        invdiff = 1.0 / ( dists(pos_idx(i)) - dists(neg_idx(1)) )
        w0 = -dists(neg_idx(1)) * invdiff
        w1 =  dists(pos_idx(i)) * invdiff
        tets_tmp_buf(n_tets_tmp)%v(:, pos_idx(i)) = &
           w0 * tets_tmp_buf(n_tets_tmp)%v(:, pos_idx(i)) + &
           w1 * tets_tmp_buf(n_tets_tmp)%v(:, neg_idx(1))
      end do
      ! The colours will have been inherited already; we just need to zero
      ! the one corresponding to the plane cut
      tets_tmp_buf(n_tets_tmp)%colours(face_no(pos_idx(1), pos_idx(2), pos_idx(3))) = 0
    case(2)
      select case(neg_cnt)
      case(2)
        ! ++--
        do i=1,pos_cnt
          invdiff = 1.0 / ( dists(pos_idx(i)) - dists(neg_idx(1)) )
          w0 = -dists(neg_idx(1)) * invdiff
          w1 =  dists(pos_idx(i)) * invdiff
          tet_tmp%v(:, i) = w0 * tet%v(:, pos_idx(i)) + w1 * tet%v(:, neg_idx(1))
        end do
        do i=1,neg_cnt
          invdiff = 1.0 / ( dists(pos_idx(i)) - dists(neg_idx(2)) )
          w0 = -dists(neg_idx(2)) * invdiff
          w1 =  dists(pos_idx(i)) * invdiff
          tet_tmp%v(:, i+2) = w0 * tet%v(:, pos_idx(i)) + w1 * tet%v(:, neg_idx(2))
        end do

        n_tets_tmp = n_tets_tmp + 1
        tets_tmp_buf(n_tets_tmp) = tet
        tets_tmp_buf(n_tets_tmp)%v(:, pos_idx(1)) = tet_tmp%v(:, 3)
        tets_tmp_buf(n_tets_tmp)%v(:, pos_idx(2)) = tet_tmp%v(:, 2)
        tets_tmp_buf(n_tets_tmp)%colours(neg_idx(1)) = 0
        tets_tmp_buf(n_tets_tmp)%colours(neg_idx(2)) = 0

        n_tets_tmp = n_tets_tmp + 1
        tets_tmp_buf(n_tets_tmp)%v(:, 1) = tet%v(:, neg_idx(2))
        tets_tmp_buf(n_tets_tmp)%colours(1) = 0
        tets_tmp_buf(n_tets_tmp)%v(:, 2) = tet_tmp%v(:, 4)
        tets_tmp_buf(n_tets_tmp)%colours(2) = 0
        tets_tmp_buf(n_tets_tmp)%v(:, 3) = tet_tmp%v(:, 3)
        tets_tmp_buf(n_tets_tmp)%colours(3) = tet%colours(pos_idx(1))
        tets_tmp_buf(n_tets_tmp)%v(:, 4) = tet_tmp%v(:, 2)
        tets_tmp_buf(n_tets_tmp)%colours(4) = tet%colours(neg_idx(1))

        n_tets_tmp = n_tets_tmp + 1
        tets_tmp_buf(n_tets_tmp)%v(:, 1) = tet%v(:, neg_idx(1))
        tets_tmp_buf(n_tets_tmp)%colours(1) = 0
        tets_tmp_buf(n_tets_tmp)%v(:, 2) = tet_tmp%v(:, 1)
        tets_tmp_buf(n_tets_tmp)%colours(2) = 0
        tets_tmp_buf(n_tets_tmp)%v(:, 3) = tet_tmp%v(:, 2)
        tets_tmp_buf(n_tets_tmp)%colours(3) = tet%colours(pos_idx(2))
        tets_tmp_buf(n_tets_tmp)%v(:, 4) = tet_tmp%v(:, 3)
        tets_tmp_buf(n_tets_tmp)%colours(4) = tet%colours(neg_idx(2))
      case(1)
        ! ++-0
        n_tets_tmp = n_tets_tmp + 1
        tets_tmp_buf(n_tets_tmp) = tet
        do i=1,pos_cnt
          invdiff = 1.0 / ( dists(pos_idx(i)) - dists(neg_idx(1)) )
          w0 = -dists(neg_idx(1)) * invdiff
          w1 =  dists(pos_idx(i)) * invdiff
          tets_tmp_buf(n_tets_tmp)%v(:, pos_idx(i)) = &
             w0 * tets_tmp_buf(n_tets_tmp)%v(:, pos_idx(i)) + &
             w1 * tets_tmp_buf(n_tets_tmp)%v(:, neg_idx(1))
        end do
        tets_tmp_buf(n_tets_tmp)%colours(neg_idx(1)) = 0
      end select
    case(1)
      select case(neg_cnt)
      case(3)
        ! +---
        do i=1,neg_cnt
          invdiff = 1.0 / ( dists(pos_idx(1)) - dists(neg_idx(i)) )
          w0 = -dists(neg_idx(i)) * invdiff
          w1 =  dists(pos_idx(1)) * invdiff
          tet_tmp%v(:, i) = w0 * tet%v(:, pos_idx(1)) + w1 * tet%v(:, neg_idx(i))
        end do

        n_tets_tmp = n_tets_tmp + 1
        tets_tmp_buf(n_tets_tmp) = tet
        tets_tmp_buf(n_tets_tmp)%v(:, pos_idx(1)) = tet_tmp%v(:, 1)
        tets_tmp_buf(n_tets_tmp)%colours(neg_idx(1)) = 0

        n_tets_tmp = n_tets_tmp + 1
        tets_tmp_buf(n_tets_tmp)%v(:, 1) = tet_tmp%v(:, 1)
        tets_tmp_buf(n_tets_tmp)%colours(1) = tet%colours(neg_idx(1))
        tets_tmp_buf(n_tets_tmp)%v(:, 2) = tet%v(:, neg_idx(2))
        tets_tmp_buf(n_tets_tmp)%colours(2) = 0
        tets_tmp_buf(n_tets_tmp)%v(:, 3) = tet%v(:, neg_idx(3))
        tets_tmp_buf(n_tets_tmp)%colours(3) = tet%colours(neg_idx(3))
        tets_tmp_buf(n_tets_tmp)%v(:, 4) = tet_tmp%v(:, 2)
        tets_tmp_buf(n_tets_tmp)%colours(4) = 0

        n_tets_tmp = n_tets_tmp + 1
        tets_tmp_buf(n_tets_tmp)%v(:, 1) = tet%v(:, neg_idx(3))
        tets_tmp_buf(n_tets_tmp)%colours(1) = 0
        tets_tmp_buf(n_tets_tmp)%v(:, 2) = tet_tmp%v(:, 2)
        tets_tmp_buf(n_tets_tmp)%colours(2) = tet%colours(neg_idx(2))
        tets_tmp_buf(n_tets_tmp)%v(:, 3) = tet_tmp%v(:, 3)
        tets_tmp_buf(n_tets_tmp)%colours(3) = 0
        tets_tmp_buf(n_tets_tmp)%v(:, 4) = tet_tmp%v(:, 1)
        tets_tmp_buf(n_tets_tmp)%colours(4) = tet%colours(neg_idx(1))
      case(2)
        ! +--0
        do i=1,neg_cnt
          invdiff = 1.0 / ( dists(pos_idx(1)) - dists(neg_idx(i)) )
          w0 = -dists(neg_idx(i)) * invdiff
          w1 =  dists(pos_idx(1)) * invdiff
          tet_tmp%v(:, i) = w0 * tet%v(:, pos_idx(1)) + w1 * tet%v(:, neg_idx(i))
        end do

        n_tets_tmp = n_tets_tmp + 1
        tets_tmp_buf(n_tets_tmp) = tet
        tets_tmp_buf(n_tets_tmp)%v(:, pos_idx(1)) = tet_tmp%v(:, 1)
        tets_tmp_buf(n_tets_tmp)%colours(neg_idx(1)) = 0

        n_tets_tmp = n_tets_tmp + 1
        tets_tmp_buf(n_tets_tmp)%v(:, 1) = tet_tmp%v(:, 2)
        tets_tmp_buf(n_tets_tmp)%colours(1) = 0
        tets_tmp_buf(n_tets_tmp)%v(:, 2) = tet%v(:, zer_idx(1))
        tets_tmp_buf(n_tets_tmp)%colours(2) = tet%colours(zer_idx(1))
        tets_tmp_buf(n_tets_tmp)%v(:, 3) = tet%v(:, neg_idx(2))
        tets_tmp_buf(n_tets_tmp)%colours(3) = 0
        tets_tmp_buf(n_tets_tmp)%v(:, 4) = tet_tmp%v(:, 1)
        tets_tmp_buf(n_tets_tmp)%colours(4) = tet%colours(neg_idx(1))
      case(1)
        ! +-00
        invdiff = 1.0 / ( dists(pos_idx(1)) - dists(neg_idx(1)) )
        w0 = -dists(neg_idx(1)) * invdiff
        w1 =  dists(pos_idx(1)) * invdiff

        n_tets_tmp = n_tets_tmp + 1
        tets_tmp_buf(n_tets_tmp) = tet
        tets_tmp_buf(n_tets_tmp)%v(:, pos_idx(1)) = w0 * tet%v(:, pos_idx(1)) + w1 * tet%v(:, neg_idx(1))
        tets_tmp_buf(n_tets_tmp)%colours(neg_idx(1)) = 0
      end select
    end select

  end subroutine clip_buf

  subroutine clip(plane, tet, tets_new, ntets_new)
    ! Clip tet against the plane.
    type(plane_type), intent(in) :: plane
    type(tet_type), intent(in) :: tet
    type(tet_type), dimension(:), intent(inout) :: tets_new
    integer, intent(out) :: ntets_new

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
      ntets_new = ntets_new + 1
      tets_new(ntets_new) = tet
      return
    end if

    ! The tet is split by the plane, so we have more work to do.

    select case(pos_cnt)
    case(3)
      ! +++-
      ntets_new = ntets_new + 1
      tets_new(ntets_new) = tet
      do i=1,pos_cnt
        invdiff = 1.0 / ( dists(pos_idx(i)) - dists(neg_idx(1)) )
        w0 = -dists(neg_idx(1)) * invdiff
        w1 =  dists(pos_idx(i)) * invdiff
        tets_new(ntets_new)%v(:, pos_idx(i)) = &
           w0 * tets_new(ntets_new)%v(:, pos_idx(i)) + &
           w1 * tets_new(ntets_new)%v(:, neg_idx(1))
      end do
      ! The colours will have been inherited already; we just need to zero
      ! the one corresponding to the plane cut
      tets_new(ntets_new)%colours(face_no(pos_idx(1), pos_idx(2), pos_idx(3))) = 0
    case(2)
      select case(neg_cnt)
      case(2)
        ! ++--
        do i=1,pos_cnt
          invdiff = 1.0 / ( dists(pos_idx(i)) - dists(neg_idx(1)) )
          w0 = -dists(neg_idx(1)) * invdiff
          w1 =  dists(pos_idx(i)) * invdiff
          tet_tmp%v(:, i) = w0 * tet%v(:, pos_idx(i)) + w1 * tet%v(:, neg_idx(1))
        end do
        do i=1,neg_cnt
          invdiff = 1.0 / ( dists(pos_idx(i)) - dists(neg_idx(2)) )
          w0 = -dists(neg_idx(2)) * invdiff
          w1 =  dists(pos_idx(i)) * invdiff
          tet_tmp%v(:, i+2) = w0 * tet%v(:, pos_idx(i)) + w1 * tet%v(:, neg_idx(2))
        end do

        ntets_new = ntets_new + 1
        tets_new(ntets_new) = tet
        tets_new(ntets_new)%v(:, pos_idx(1)) = tet_tmp%v(:, 3)
        tets_new(ntets_new)%v(:, pos_idx(2)) = tet_tmp%v(:, 2)
        tets_new(ntets_new)%colours(neg_idx(1)) = 0
        tets_new(ntets_new)%colours(neg_idx(2)) = 0

        ntets_new = ntets_new + 1
        tets_new(ntets_new)%v(:, 1) = tet%v(:, neg_idx(2))
        tets_new(ntets_new)%colours(1) = 0
        tets_new(ntets_new)%v(:, 2) = tet_tmp%v(:, 4)
        tets_new(ntets_new)%colours(2) = 0
        tets_new(ntets_new)%v(:, 3) = tet_tmp%v(:, 3)
        tets_new(ntets_new)%colours(3) = tet%colours(pos_idx(1))
        tets_new(ntets_new)%v(:, 4) = tet_tmp%v(:, 2)
        tets_new(ntets_new)%colours(4) = tet%colours(neg_idx(1))

        ntets_new = ntets_new + 1
        tets_new(ntets_new)%v(:, 1) = tet%v(:, neg_idx(1))
        tets_new(ntets_new)%colours(1) = 0
        tets_new(ntets_new)%v(:, 2) = tet_tmp%v(:, 1)
        tets_new(ntets_new)%colours(2) = 0
        tets_new(ntets_new)%v(:, 3) = tet_tmp%v(:, 2)
        tets_new(ntets_new)%colours(3) = tet%colours(pos_idx(2))
        tets_new(ntets_new)%v(:, 4) = tet_tmp%v(:, 3)
        tets_new(ntets_new)%colours(4) = tet%colours(neg_idx(2))
      case(1)
        ! ++-0
        ntets_new = ntets_new + 1
        tets_new(ntets_new) = tet
        do i=1,pos_cnt
          invdiff = 1.0 / ( dists(pos_idx(i)) - dists(neg_idx(1)) )
          w0 = -dists(neg_idx(1)) * invdiff
          w1 =  dists(pos_idx(i)) * invdiff
          tets_new(ntets_new)%v(:, pos_idx(i)) = &
             w0 * tets_new(ntets_new)%v(:, pos_idx(i)) + &
             w1 * tets_new(ntets_new)%v(:, neg_idx(1))
        end do
        tets_new(ntets_new)%colours(neg_idx(1)) = 0
      end select
    case(1)
      select case(neg_cnt)
      case(3)
        ! +---
        do i=1,neg_cnt
          invdiff = 1.0 / ( dists(pos_idx(1)) - dists(neg_idx(i)) )
          w0 = -dists(neg_idx(i)) * invdiff
          w1 =  dists(pos_idx(1)) * invdiff
          tet_tmp%v(:, i) = w0 * tet%v(:, pos_idx(1)) + w1 * tet%v(:, neg_idx(i))
        end do

        ntets_new = ntets_new + 1
        tets_new(ntets_new) = tet
        tets_new(ntets_new)%v(:, pos_idx(1)) = tet_tmp%v(:, 1)
        tets_new(ntets_new)%colours(neg_idx(1)) = 0

        ntets_new = ntets_new + 1
        tets_new(ntets_new)%v(:, 1) = tet_tmp%v(:, 1)
        tets_new(ntets_new)%colours(1) = tet%colours(neg_idx(1))
        tets_new(ntets_new)%v(:, 2) = tet%v(:, neg_idx(2))
        tets_new(ntets_new)%colours(2) = 0
        tets_new(ntets_new)%v(:, 3) = tet%v(:, neg_idx(3))
        tets_new(ntets_new)%colours(3) = tet%colours(neg_idx(3))
        tets_new(ntets_new)%v(:, 4) = tet_tmp%v(:, 2)
        tets_new(ntets_new)%colours(4) = 0

        ntets_new = ntets_new + 1
        tets_new(ntets_new)%v(:, 1) = tet%v(:, neg_idx(3))
        tets_new(ntets_new)%colours(1) = 0
        tets_new(ntets_new)%v(:, 2) = tet_tmp%v(:, 2)
        tets_new(ntets_new)%colours(2) = tet%colours(neg_idx(2))
        tets_new(ntets_new)%v(:, 3) = tet_tmp%v(:, 3)
        tets_new(ntets_new)%colours(3) = 0
        tets_new(ntets_new)%v(:, 4) = tet_tmp%v(:, 1)
        tets_new(ntets_new)%colours(4) = tet%colours(neg_idx(1))
      case(2)
        ! +--0
        do i=1,neg_cnt
          invdiff = 1.0 / ( dists(pos_idx(1)) - dists(neg_idx(i)) )
          w0 = -dists(neg_idx(i)) * invdiff
          w1 =  dists(pos_idx(1)) * invdiff
          tet_tmp%v(:, i) = w0 * tet%v(:, pos_idx(1)) + w1 * tet%v(:, neg_idx(i))
        end do

        ntets_new = ntets_new + 1
        tets_new(ntets_new) = tet
        tets_new(ntets_new)%v(:, pos_idx(1)) = tet_tmp%v(:, 1)
        tets_new(ntets_new)%colours(neg_idx(1)) = 0

        ntets_new = ntets_new + 1
        tets_new(ntets_new)%v(:, 1) = tet_tmp%v(:, 2)
        tets_new(ntets_new)%colours(1) = 0
        tets_new(ntets_new)%v(:, 2) = tet%v(:, zer_idx(1))
        tets_new(ntets_new)%colours(2) = tet%colours(zer_idx(1))
        tets_new(ntets_new)%v(:, 3) = tet%v(:, neg_idx(2))
        tets_new(ntets_new)%colours(3) = 0
        tets_new(ntets_new)%v(:, 4) = tet_tmp%v(:, 1)
        tets_new(ntets_new)%colours(4) = tet%colours(neg_idx(1))
      case(1)
        ! +-00
        invdiff = 1.0 / ( dists(pos_idx(1)) - dists(neg_idx(1)) )
        w0 = -dists(neg_idx(1)) * invdiff
        w1 =  dists(pos_idx(1)) * invdiff

        ntets_new = ntets_new + 1
        tets_new(ntets_new) = tet
        tets_new(ntets_new)%v(:, pos_idx(1)) = w0 * tet%v(:, pos_idx(1)) + w1 * tet%v(:, neg_idx(1))
        tets_new(ntets_new)%colours(neg_idx(1)) = 0
      end select
    end select

  end subroutine clip
  
  pure function get_planes_tet(tet) result(plane)
    type(tet_type), intent(in) :: tet
    type(plane_type), dimension(4) :: plane

    real, dimension(3) :: edge10, edge20, edge30, edge21, edge31
    real :: det
    integer :: i

    edge10 = tet%v(:, 2) - tet%v(:, 1);
    edge20 = tet%v(:, 3) - tet%v(:, 1);
    edge30 = tet%v(:, 4) - tet%v(:, 1);
    edge21 = tet%v(:, 3) - tet%v(:, 2);
    edge31 = tet%v(:, 4) - tet%v(:, 2);

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
      plane(i)%c = dot_product(tet%v(:, i), plane(i)%normal)
    end do

  end function get_planes_tet

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
      dists(i) = dot_product(plane%normal, tet%v(:, i)) - plane%c
    end forall
  end function distances_to_plane

  pure function tetrahedron_volume_real(tet) result(volume)
    real, dimension(3, 4), intent(in) :: tet

    real :: volume

    real, dimension(3) :: e1, e2, e3

    e1 = tet(:, 2) - tet(:, 1)
    e2 = tet(:, 3) - tet(:, 1)
    e3 = tet(:, 4) - tet(:, 1)

    volume = (1.0 / 6.0) * abs(e1(1) * (e2(2) * e3(3) - e2(3) * e3(2)) &
                           & + e1(2) * (e2(3) * e3(1) - e2(1) * e3(3)) &
                           & + e1(3) * (e2(1) * e3(2) - e2(2) * e3(1)))

  end function tetrahedron_volume_real

  pure function tetrahedron_volume_tet(tet) result(vol)
    type(tet_type), intent(in) :: tet
    real :: vol
    real, dimension(3) :: cross, vecA, vecB, vecC

    vecA = tet%v(:, 1) - tet%v(:, 4)
    vecB = tet%v(:, 2) - tet%v(:, 4)
    vecC = tet%v(:, 3) - tet%v(:, 4)

    cross(1) = vecB(2) * vecC(3) - vecB(3) * vecC(2)
    cross(2) = vecB(3) * vecC(1) - vecB(1) * vecC(3)
    cross(3) = vecB(1) * vecC(2) - vecB(2) * vecC(1)

    vol = abs(dot_product(vecA, cross)) / 6.0
    
  end function tetrahedron_volume_tet

  pure function face_no(i, j, k) result(face)
    ! Given three local node numbers, what is the face that they share?
    integer, intent(in) :: i, j, k
    integer :: face

    do face=1,4
      if (face /= i .and. face /= j .and. face /= k) return
    end do

  end function face_no

end module libsupermesh_tet_intersection
