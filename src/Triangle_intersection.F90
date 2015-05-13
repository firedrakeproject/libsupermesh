#define BUF_SIZE 150
#include "fdebug.h"

!
! Original Triangle (A)
!
!         E(Ex, Ey)
!         /\
!        /  \
!       /    \
!      /      \
!     /--------\
!  N(Nx,Ny)     D(Dx,Dy)
!
! <==================================>
!
! Clipping triangle (B)
!
!       Q(Qx, Qy)
!       /\
!      /  \
!     /    \
!    /      \
!   /--------\
! P(Px,Py)    R(Rx, Ry)
!
! <==================================>
!
! Intersecting Triangles
!
!      Q(Qx,Qy)  E(Ex, Ey)
!            |   /\
!  <---------|  /  \
!  ~         | /    \
!  n         |/      \
!            *        \
!           /|         \
!          / |          \
!         /  |           \
!        /   |            \
!       /    |             \
!      /     |              \
!     /      |I(Ix,Iy)       \
!    /-------*----------------\
!   N(Nx,Ny) |            D(Dx,Dy)
!            |
!            |
!O(Ox,Oy)    P(Px,Py)
!
! Steps:
! a. Find the normals of all sides of Triangle B.
! b. Transform the normals so that their direction is consistent.
!                                                    ->   ~
! c. Decide if we want to keep a node of Triangle A (PN * n >0).
!                                ->    ->   ->       ->         ->   ->
! d. Find the co-ordinates of I (NI = -ON + OI = l * ND = l * (-ON + OD)


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
    real, dimension(2, 3) :: V ! vertices of the tri
    integer, dimension(3) :: colours = -1 ! surface colours
  end type tri_type

  type line_type
    real, dimension(2) :: normal
    real :: c = 0
  end type line_type

  type(tri_type), dimension(BUF_SIZE), save :: tri_array, tri_array_tmp
  integer :: tri_cnt = 0, tri_cnt_tmp = 0
  type(mesh_type), save :: intersection_mesh
  logical, save :: mesh_allocated = .false.

  public :: tri_type, line_type, libsupermesh_intersect_tris, get_lines, finalise_tri_intersector

  interface libsupermesh_intersect_tris
    module procedure libsupermesh_intersect_tris_dt, libsupermesh_intersect_tris_dt_public
  end interface

  interface get_lines
    module procedure get_lines_tri
  end interface

  contains

  subroutine finalise_tri_intersector
    if (mesh_allocated) then
      call deallocate(intersection_mesh)
      mesh_allocated = .false.
    end if
  end subroutine finalise_tri_intersector

  subroutine libsupermesh_intersect_tris_dt_public(triA_V, triA_colours, &
     sizeOfLinesB, linesB_normal, linesB_c, &
     quadVertices, quadDim, quadNgi, quadDegree, shapeLoc, shapeDim, shapeDegree, &
     stat, output, &
     surface_shape, surface_positions, surface_colours)
    real, intent(in), dimension(2, 3) :: triA_V
    integer, intent(in), dimension(3) :: triA_colours
    integer, intent(in)               :: sizeOfLinesB
    real, intent(in), dimension(sizeOfLinesB,2) :: linesB_normal
    real, intent(in), dimension(sizeOfLinesB)   :: linesB_c
    type(vector_field), intent(inout) :: output
    type(vector_field), intent(out), optional :: surface_positions
    type(scalar_field), intent(out), optional :: surface_colours
    type(element_type), intent(in), optional :: surface_shape
    integer, intent(in) :: quadVertices, quadDim, quadNgi, quadDegree, shapeLoc, shapeDim, shapeDegree
    integer :: ele, i, j, k, l
    integer, intent(out) :: stat

    type(tri_type)  :: triA
    type(line_type), dimension(sizeOfLinesB) :: linesB
    
    if (present(surface_colours) .or. present(surface_positions) .or. present(surface_shape)) then
      assert(present(surface_positions))
      assert(present(surface_colours))
      assert(present(surface_shape))
    end if

    triA%V = triA_V
    triA%colours = triA_colours

    do i = 1, size(linesB)
      do k = 1, 3
        linesB(i)%normal(k)=linesB_normal(i,k)
      end do
    end do
    FORALL(i=1:sizeOfLinesB) linesB(i)%c=linesB_c(i)
    
    assert(shapeDegree == 1)
    assert(shapeNumberingFamily == FAMILY_SIMPLEX)
    assert(shapeDim == 3)
    
    output = libsupermesh_intersect_tris_dt(triA, linesB, &
        quadVertices, quadDim, quadNgi, quadDegree, &
        shapeLoc, shapeDim, shapeDegree, &
        stat = stat, output = output)
  
  end subroutine libsupermesh_intersect_tris_dt_public
  
!  subroutine libsupermesh_intersect_tris_dt(triA, linesB, shape, stat, output, surface_shape, surface_positions, surface_colours)
  subroutine libsupermesh_intersect_tris_dt(triA, linesB, &
    quadVertices, quadDim, quadNgi, quadDegree, shapeLoc, shapeDim, shapeDegree, &
    stat, output, surface_shape, surface_positions, surface_colours)
    type(tri_type), intent(in) :: triA
    type(line_type), dimension(:), intent(in) :: linesB
    type(vector_field), intent(inout) :: output
    type(vector_field), intent(out), optional :: surface_positions
    type(scalar_field), intent(out), optional :: surface_colours
    type(element_type), intent(in), optional :: surface_shape
    integer, intent(in) :: quadVertices, quadDim, quadNgi, quadDegree, shapeLoc, shapeDim, shapeDegree
    integer :: ele, i, j, k, l
    integer, intent(out) :: stat

!    real, dimension(3) :: vec_tmp
!    integer, dimension(3) :: idx_tmp
!    integer :: surface_eles, colour_tmp
!    type(mesh_type) :: surface_mesh, pwc_surface_mesh

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
!      call allocate(intersection_mesh, BUF_SIZE * 3, BUF_SIZE, shape, name="IntersectionMesh")
      call allocate(intersection_mesh, BUF_SIZE * 3, BUF_SIZE, shapeLoc, shapeDim, shapeDegree, quadVertices, quadDim, quadNgi, quadDegree, name="IntersectionMesh")
      intersection_mesh%ndglno = (/ (i, i=1,BUF_SIZE*3) /)
      intersection_mesh%continuity = -1
      mesh_allocated = .true.
    end if

    do i=1,size(linesB)
      ! Clip the tri_array against the i'th plane
      tri_cnt_tmp = 0

      do j=1,tri_cnt
        call clip(linesB(i), tri_array(j))
      end do

      if (i /= size(linesB)) then
        tri_cnt = tri_cnt_tmp
        tri_array(1:tri_cnt) = tri_array_tmp(1:tri_cnt)
      else
        ! Copy the result
        tri_cnt = 0
        do j=1,tri_cnt_tmp
          tri_cnt = tri_cnt + 1
          tri_array(tri_cnt) = tri_array_tmp(j)
        end do
      end if
    end do

    if (tri_cnt == 0) then
      stat=1
      return
    end if

    stat = 0
    intersection_mesh%nodes = tri_cnt*3
    intersection_mesh%elements = tri_cnt
    call allocate(output, 2, intersection_mesh, "IntersectionCoordinates")

    do ele=1,tri_cnt
      call set(output, ele_nodes(output, ele), tri_array(ele)%V)
    end do

    if (present(surface_positions)) then
      FLAbort("libsupermesh_intersect_tris_dt: surface_positions are not supported.")
    end if

  end subroutine libsupermesh_intersect_tris_dt

  subroutine clip(line, tri)
  ! Clip tri against the plane
  ! and append any output to tri_array_tmp.
    type(line_type), intent(in) :: line
    type(tri_type), intent(in) :: tri

    real, dimension(3) :: dists
    integer :: neg_cnt, pos_cnt, zer_cnt
    integer, dimension(3) :: neg_idx, pos_idx, zer_idx
    type(line_type), dimension(2) :: temp_lines
    integer :: i

    real :: invdiff, w0, w1
    type(tri_type) :: tri_tmp

    ! Negative == inside
    ! Positive == outside

    neg_cnt = 0
    pos_cnt = 0
    
    dists = distances_to_line(line, tri)
    do i=1,3
      if (dists(i) < 0.0) then
        neg_cnt = neg_cnt + 1
        neg_idx(neg_cnt) = i
      else if (dists(i) > 0.0) then
        pos_cnt = pos_cnt + 1
        pos_idx(pos_cnt) = i
      end if
    end do

    if (pos_cnt == 0) then
      ! tri is completely on positive side of line, full clip
      return
    end if

    if (neg_cnt == 0) then
      ! tri is completely on negative side of line, no clip
      tri_cnt_tmp = tri_cnt_tmp + 1
      tri_array_tmp(tri_cnt_tmp) = tri
      return
    end if

    ! The tri is split by the line, so we have more work to do.

    select case(pos_cnt)
    case(1)
      select case(neg_cnt)
      case(2)
        tri_cnt_tmp = tri_cnt_tmp + 1
        tri_array_tmp(tri_cnt_tmp) = tri

        do i=1,neg_cnt
          invdiff = 1.0 / ( dists(pos_idx(1)) - dists(neg_idx(i)) )
          w0 = -dists(neg_idx(i)) * invdiff
          w1 =  dists(pos_idx(1)) * invdiff
          tri_array_tmp(tri_cnt_tmp)%V(:, neg_idx(i)) = &
               w0 * tri%V(:, pos_idx(1)) + w1 * tri%V(:, neg_idx(i))
        end do
      case default
        FLAbort("Error. Found more than three points.")
      end select
    case(2)
      select case(neg_cnt)
      case(1)
        ! We need to return two triangles back
        do i=1,pos_cnt
          invdiff = 1.0 / ( dists(pos_idx(i)) - dists(neg_idx(1)) )
          w0 = -dists(neg_idx(1)) * invdiff
          w1 =  dists(pos_idx(i)) * invdiff
          tri_tmp%V(:,i) = &
               w0 * tri%V(:, pos_idx(i)) + w1 * tri%V(:, neg_idx(1))
        end do

        tri_cnt_tmp = tri_cnt_tmp + 1
        tri_array_tmp(tri_cnt_tmp) = tri
        tri_array_tmp(tri_cnt_tmp)%V(:, pos_idx(1)) = tri_tmp%V(:, 1)
        tri_array_tmp(tri_cnt_tmp)%V(:, pos_idx(2)) = tri_tmp%V(:, 2)
        tri_array_tmp(tri_cnt_tmp)%V(:, neg_idx(1)) = tri%V(:, pos_idx(1))

        tri_cnt_tmp = tri_cnt_tmp + 1
        tri_array_tmp(tri_cnt_tmp) = tri
        tri_array_tmp(tri_cnt_tmp)%V(:, pos_idx(1)) = tri%V(:, pos_idx(1))
        tri_array_tmp(tri_cnt_tmp)%V(:, pos_idx(2)) = tri_tmp%V(:, 2)
        tri_array_tmp(tri_cnt_tmp)%V(:, neg_idx(1)) = tri%V(:, pos_idx(2))
      case default
        FLAbort("Error. Found more than three points.")
      end select
    end select
  end subroutine clip

! Find the normals of all sides of triangle B (steps a and b).
!  pure function get_lines_tri(tri) result(lines)
  function get_lines_tri(tri) result(lines)
    type(tri_type), intent(in) :: tri
    type(line_type), dimension(3) :: lines

    real, dimension(2) :: edge10, edge20, edge30
    real :: det
    integer :: i
    
    ! Find the normals of all sides of triangle B (step a).
    ! P = 1, Q = 2, R = 3
    edge10 = tri%V(:, 2) - tri%V(:, 1); ! PQ = ( Qx - Px , Qy - Py )
    edge20 = tri%V(:, 3) - tri%V(:, 2); ! QR = ( Rx - Qx , Ry - Qy )
    edge30 = tri%V(:, 1) - tri%V(:, 3); ! RP = ( Px - Rx , Py - Ry )

    lines(1)%normal = unit_cross(edge10) ! PQ normal = ( -PQy , PQx )
    lines(2)%normal = unit_cross(edge20) ! QR normal = ( -QRy , QRx )
    lines(3)%normal = unit_cross(edge30) ! PR normal = ( -PRy , PRx )

    ! Transform the normals so that their direction is consistent (step b).
    det = dot_product(edge10, lines(3)%normal)
    if (det < 0) then
      do i=1,3
        lines(i)%normal = -lines(i)%normal
      end do
    end if

    ! And calibrate what is the zero of this line by dotting with
    ! a point we know to be on it
    do i=1,3
      lines(i)%c = dot_product(tri%V(:, i), lines(i)%normal)
    end do

  end function get_lines_tri

! Return a cross product (in 2D).  
  pure function unit_cross(vecA) result(cross)
    real, dimension(2), intent(in) :: vecA
    real, dimension(2) :: cross
    cross(1) = (-1) * vecA(2)
    cross(2) = vecA(1)

  end function unit_cross

  pure function distances_to_line(line, tri) result(dists)
    type(line_type), intent(in) :: line
    type(tri_type), intent(in) :: tri
    real, dimension(3) :: dists
    real, dimension(2) :: P
    integer :: i

    forall(i=1:3)
      dists(i) = dot_product(line%normal, tri%V(:, i)) - line%c
    end forall

  end function distances_to_line

end module libsupermesh_tri_intersection_module
