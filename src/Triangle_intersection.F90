#define BUF_SIZE 150
#include "fdebug.h"

module libsupermesh_tri_intersection_module

  use libsupermesh_elements
  use libsupermesh_vector_tools
  use libsupermesh_fields_data_types
  use libsupermesh_fields_base
  use libsupermesh_fields_allocates
  use libsupermesh_fields_manipulation
  use libsupermesh_transform_elements
  implicit none

  type tri_type
    real, dimension(3, 3) :: V ! vertices of the tri
    integer, dimension(3) :: colours = -1 ! surface colours
  end type tri_type

  type plane_type
    real, dimension(3) :: normal
    real :: c
  end type plane_type

  type(tri_type), dimension(BUF_SIZE), save :: tri_array, tri_array_tmp
  integer :: tri_cnt = 0, tri_cnt_tmp = 0
  type(mesh_type), save :: intersection_mesh
  logical, save :: mesh_allocated = .false.

  public :: tri_type, plane_type, libsupermesh_intersect_tris, get_planes, finalise_tri_intersector

  interface libsupermesh_intersect_tris
    module procedure libsupermesh_intersect_tris_dt
  end interface

  interface get_planes
    module procedure get_planes_tri
  end interface

  contains

  subroutine finalise_tri_intersector
    if (mesh_allocated) then
      call deallocate(intersection_mesh)
      mesh_allocated = .false.
    end if
  end subroutine finalise_tri_intersector

  subroutine libsupermesh_intersect_tris_dt(triA, planesB, shape, stat, output, surface_shape, surface_positions, surface_colours)
    type(tri_type), intent(in) :: triA
    type(plane_type), dimension(:), intent(in) :: planesB
    type(element_type), intent(in) :: shape
    type(vector_field), intent(inout) :: output
    type(vector_field), intent(out), optional :: surface_positions
    type(scalar_field), intent(out), optional :: surface_colours
    type(element_type), intent(in), optional :: surface_shape
    integer :: ele
    integer, intent(out) :: stat

    integer :: i, j, k, l
    real :: vol
    real, dimension(3) :: vec_tmp
    integer, dimension(3) :: idx_tmp
    integer :: surface_eles, colour_tmp
    type(mesh_type) :: surface_mesh, pwc_surface_mesh

    if (present(surface_colours) .or. present(surface_positions) .or. present(surface_shape)) then
      assert(present(surface_positions))
      assert(present(surface_colours))
      assert(present(surface_shape))
    end if


    assert(shape%degree == 1)
    assert(shape%numbering%family == FAMILY_SIMPLEX)
    assert(shape%dim == 3)
    
    tri_cnt = 1
    tri_array(1) = triA

    if (.not. mesh_allocated) then
      call allocate(intersection_mesh, BUF_SIZE * 4, BUF_SIZE, shape, name="IntersectionMesh")
      intersection_mesh%ndglno = (/ (i, i=1,BUF_SIZE*4) /)
      intersection_mesh%continuity = -1
      mesh_allocated = .true.
    end if
    
    do i=1,size(planesB)
      ! Clip the tri_array against the i'th plane
!      tri_cnt_tmp = 0

      do j=1,tri_cnt
        call clip(planesB(i), tri_array(j))
        tri_cnt = tri_cnt + 1
        tri_array(tri_cnt) = tri_array_tmp(j)
      end do
    end do
    
    stat = 0
    intersection_mesh%nodes = tri_cnt*4
    intersection_mesh%elements = tri_cnt
    call allocate(output, 3, intersection_mesh, "IntersectionCoordinates")
    
    do ele=1,tri_cnt
      call set(output, ele_nodes(output, ele), tri_array(ele)%V)
    end do
    
    if (present(surface_positions)) then
      FLAbort("libsupermesh_intersect_tris_dt: Not surface_positions")
    end if

  end subroutine libsupermesh_intersect_tris_dt


  subroutine clip(plane, tri)
  ! Clip tri against the plane
  ! and append any output to tet_array_tmp.
    type(plane_type), intent(in) :: plane
    type(tri_type), intent(in) :: tri

    real, dimension(3) :: dists
    integer :: neg_cnt, pos_cnt, zer_cnt
    integer, dimension(3) :: neg_idx, pos_idx, zer_idx
    integer :: i

    real :: invdiff, w0, w1
    type(tri_type) :: tet_tmp

    ! Negative == inside
    ! Positive == outside

    neg_cnt = 0
    pos_cnt = 0
    zer_cnt = 0

  end subroutine clip

  pure function get_planes_tri(tri) result(plane)
    type(tri_type), intent(in) :: tri
    type(plane_type), dimension(4) :: plane

    real, dimension(3) :: edge10, edge20, edge30, edge21, edge31
    real :: det
    integer :: i

  end function get_planes_tri

  pure function unit_cross(vecA, vecB) result(cross)
    real, dimension(3), intent(in) :: vecA, vecB
    real, dimension(3) :: cross
    cross(1) = vecA(2) * vecB(3) - vecA(3) * vecB(2)
    cross(2) = vecA(3) * vecB(1) - vecA(1) * vecB(3)
    cross(3) = vecA(1) * vecB(2) - vecA(2) * vecB(1)

    cross = cross / norm2(cross)
  end function unit_cross

  pure function distances_to_plane(plane, tri) result(dists)
    type(plane_type), intent(in) :: plane
    type(tri_type), intent(in) :: tri
    real, dimension(3) :: dists
    integer :: i

  end function distances_to_plane

  pure function tri_volume(tet) result(vol)
    type(tri_type), intent(in) :: tet
    real :: vol
    real, dimension(3) :: cross, vecA, vecB, vecC

  end function tri_volume

  function face_no(i, j, k) result(face)
    ! Given three local node numbers, what is the face that they share?
    integer, intent(in) :: i, j, k
    integer :: face

    do face=1,4
      if (face /= i .and. face /= j .and. face /= k) return
    end do

  end function face_no

end module libsupermesh_tri_intersection_module
