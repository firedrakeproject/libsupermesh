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
    real, dimension(2, 3) :: V ! vertices of the tri
    integer, dimension(3) :: colours = -1 ! surface colours
  end type tri_type

  type line_type
    real, dimension(2) :: normal
!                 (A,B)
    real :: c
  end type line_type

  type(tri_type), dimension(BUF_SIZE), save :: tri_array, tri_array_tmp
  integer :: tri_cnt = 0, tri_cnt_tmp = 0
  type(mesh_type), save :: intersection_mesh
  logical, save :: mesh_allocated = .false.

  public :: tri_type, line_type, libsupermesh_intersect_tris, get_lines, finalise_tri_intersector

  interface libsupermesh_intersect_tris
    module procedure libsupermesh_intersect_tris_dt
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

  subroutine libsupermesh_intersect_tris_dt(triA, linesB, shape, stat, output, surface_shape, surface_positions, surface_colours)
    type(tri_type), intent(in) :: triA
    type(line_type), dimension(:), intent(in) :: linesB
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
    type(line_type), dimension(3) :: linesA

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
    linesA = get_lines(triA)

    if (.not. mesh_allocated) then
      call allocate(intersection_mesh, BUF_SIZE * 3, BUF_SIZE, shape, name="IntersectionMesh")
      intersection_mesh%ndglno = (/ (i, i=1,BUF_SIZE*3) /)
      intersection_mesh%continuity = -1
      mesh_allocated = .true.
    end if
    
    ! IAKOVOS REMOVE COMMENT
    write(*,*) "libsupermesh_intersect_tris_dt: 2: size(linesB):",size(linesB),", tri_cnt:",tri_cnt,"."
    
    do i=1,size(linesB)
      ! Clip the tri_array against the i'th plane
      tri_cnt_tmp = 0
      ! IAKOVOS REMOVE COMMENT
      write(*,*) "libsupermesh_intersect_tris_dt: 3: tri_cnt_tmp:",tri_cnt_tmp,", tri_cnt:",tri_cnt,"."

      do j=1,tri_cnt
        call clip(linesB(i), tri_array(j), linesA)
      end do
      
      if (i /= size(linesB)) then
        ! IAKOVOS REMOVE COMMENT
        write(*,*) "libsupermesh_intersect_tris_dt: Inside if i:",i,"."
        tri_cnt = tri_cnt_tmp
        tri_array(1:tri_cnt) = tri_array_tmp(1:tri_cnt)
      end if
      
!      tri_cnt = tri_cnt + 1
!      tri_array(tri_cnt) = tri_array_tmp(j)
    end do
    
    ! IAKOVOS REMOVE COMMENT
    write(*,*) "libsupermesh_intersect_tris_dt: 4"
    
    if (tri_cnt == 0) then
      stat=1
      return
    end if
    
    ! IAKOVOS REMOVE COMMENT
    write(*,*) "libsupermesh_intersect_tris_dt: 5"
    
    stat = 0
    intersection_mesh%nodes = tri_cnt*4
    intersection_mesh%elements = tri_cnt
    call allocate(output, 2, intersection_mesh, "IntersectionCoordinates")
    
    ! IAKOVOS REMOVE COMMENT
    write(*,*) "libsupermesh_intersect_tris_dt: 6"
    
    do ele=1,tri_cnt
      call set(output, ele_nodes(output, ele), tri_array(ele)%V)
    end do
    
    ! IAKOVOS REMOVE COMMENT
    write(*,*) "libsupermesh_intersect_tris_dt: 7"
    
    if (present(surface_positions)) then
      FLAbort("libsupermesh_intersect_tris_dt: Not surface_positions")
    end if

  end subroutine libsupermesh_intersect_tris_dt


  subroutine clip(line, tri, linesA)
  ! Clip tri against the plane
  ! and append any output to tri_array_tmp.
    type(line_type), intent(in) :: line
    type(tri_type), intent(in) :: tri
    type(line_type), dimension(:), intent(in) :: linesA

    real, dimension(3) :: dists
    integer :: neg_cnt, pos_cnt, zer_cnt
    integer, dimension(3) :: neg_idx, pos_idx, zer_idx
    type(line_type), dimension(2) :: temp_lines
    integer :: i

    real :: invdiff, w0, w1
    type(tri_type) :: tet_tmp

    ! Negative == ouside
    ! Positive == inside

    neg_cnt = 0
    pos_cnt = 0
    zer_cnt = 0
    
    ! IAKOVOS REMOVE COMMENT
    write(*,*) "clip: line: normal:",line%normal,", line%c:",line%c,"."
    write(*,*) "clip:         tri%V:",tri%V,"."

    dists = distances_to_line(line, tri)
    ! IAKOVOS REMOVE COMMENT
    write(*,*) "clip: distances_to_line: dists(1):",dists(1),", dists(2):",dists(2),", dists(3):",dists(3),"."
    do i=1,3
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
    
    if (pos_cnt == 0) then
      ! tri is completely on positive side of line, full clip
      ! IAKOVOS REMOVE COMMENT
      write(*,*) "clip: tri is completely on positive side of line, full clip"
      return
    end if

    if (neg_cnt == 0) then
      ! tri is completely on negative side of line, no clip
      ! IAKOVOS REMOVE COMMENT
      write(*,*) "clip: tri is completely on negative side of line, no clip"
      tri_cnt_tmp = tri_cnt_tmp + 1
      tri_array_tmp(tri_cnt_tmp) = tri
      return
    end if
    
    ! The tri is split by the line, so we have more work to do.
    
    select case(pos_cnt)
    case(1)
      ! IAKOVOS REMOVE COMMENT
      write(*,*) "clip: Case 1"
      select case(neg_cnt)
      case(2)
        ! Easy, we just need to return one triangle back.
        
        ! Find the points where the line pos_idx(pos_cnt) intersects the other two lines and return that triangle
        write(*,*) "clip: Case 1: subcase 2:"
        tri_cnt_tmp = tri_cnt_tmp + 1
        tri_array_tmp(tri_cnt_tmp) = tri
        do i=1,pos_cnt
          invdiff = 1.0 / ( dists(pos_idx(i)) - dists(neg_idx(1)) )
          write(*,*) "clip: Case 2: subcase 1:",i," point:",pos_idx(i),"(",tri%V(:,pos_idx(i)),") is POSITIVE (invdiff:",invdiff,")."
        
        end do
        do i=1,neg_cnt
          invdiff = 1.0 / ( dists(pos_idx(1)) - dists(neg_idx(i)) )
          write(*,*) "clip: Case 1: subcase 2:",i," point:",neg_idx(i)," is NEGATIVE (invdiff:",invdiff,")."
          
        end do
      case default
        FLAbort("Error. Found more than three points.")
      end select
    case(2)
      ! IAKOVOS REMOVE COMMENT
      write(*,*) "clip: Case 2"
      select case(neg_cnt)
      case(1)
        ! We need to return two triangles
        
        ! Find the points w
        write(*,*) "clip: Case 2: subcase 1:"
        
!        temp_lines
        
        tri_cnt_tmp = tri_cnt_tmp + 1
        tri_array_tmp(tri_cnt_tmp) = tri
        do i=1,pos_cnt
          invdiff = 1.0 / ( dists(pos_idx(i)) - dists(neg_idx(1)) )
          write(*,*) "clip: Case 2: subcase 1:",i," point:",pos_idx(i),"(",tri%V(:,pos_idx(i)),") is POSITIVE (invdiff:",invdiff,")."
!          if ( )
        
        end do
        do i=1,neg_cnt
          invdiff = 1.0 / ( dists(pos_idx(1)) - dists(neg_idx(i)) )
          write(*,*) "clip: Case 2: subcase 1:",i," point:",neg_idx(i)," is NEGATIVE (invdiff:",invdiff,")."
          
        end do
        
        
        
      case default
        FLAbort("Error. Found more than three points.")
      end select
    end select
    
    
    ! IAKOVOS REMOVE COMMENT
    write(*,*) "clip: ..."
    
  end subroutine clip

!  pure function get_lines_tri(tri) result(lines)
  function get_lines_tri(tri) result(lines)
    type(tri_type), intent(in) :: tri
    type(line_type), dimension(3) :: lines

    real, dimension(2) :: edge12, edge13, edge32
    real :: det
    integer :: i
    
!    do i=1,3
    ! (1,) == X; (2,) == Y
    ! 1->2
    lines(1)%normal(1) = tri%V(2,1) - tri%V(2,2) ! A = Y1 - Y2
    lines(1)%normal(2) = tri%V(1,2) - tri%V(1,1) ! B = X2 - X1
    ! 2->3
    lines(2)%normal(1) = tri%V(2,2) - tri%V(2,3) ! A = Y1 - Y2
    lines(2)%normal(2) = tri%V(1,3) - tri%V(1,2) ! B = X2 - X1
    ! 3->1
    lines(3)%normal(1) = tri%V(2,3) - tri%V(2,1) ! A = Y1 - Y2
    lines(3)%normal(2) = tri%V(1,1) - tri%V(1,3) ! B = X2 - X1
!    end do
    
!    edge12 = tri%V(:, 2) - tri%V(:, 1);
!    edge13 = tri%V(:, 3) - tri%V(:, 1);
!    edge32 = tri%V(:, 3) - tri%V(:, 2);
    
    ! IAKOVOS REMOVE COMMENT
    write(*,*) "get_lines_tri: input: point1:",tri%V(1,1),tri%V(2,1),", point2:",tri%V(1,2),tri%V(2,2),&
            ", point3:",tri%V(1,3),tri%V(2,3),"."
!    write(*,*) "get_lines_tri: edge12(1):",edge12(1),", edge12(2):",edge12(2),"."
!    write(*,*) "get_lines_tri: edge13(1):",edge13(1),", edge13(2):",edge13(2),"."
!    write(*,*) "get_lines_tri: edge32(1):",edge32(1),", edge32(2):",edge32(2),"."
 
!    lines(1)%normal = unit_cross(edge12)
!    lines(2)%normal = unit_cross(edge13)
!    lines(3)%normal = unit_cross(edge32)
    
    ! IAKOVOS REMOVE COMMENT
!    write(*,*) "get_lines_tri: lines(1)%normal:",lines(1)%normal,"."
!    write(*,*) "get_lines_tri: lines(2)%normal:",lines(2)%normal,"."
!    write(*,*) "get_lines_tri: lines(3)%normal:",lines(3)%normal,"."

!    det = dot_product(edge12, lines(3)%normal)
!    if (det < 0) then
!      do i=1,3
!        lines(i)%normal = -lines(i)%normal
!      end do
!    end if
    
    ! IAKOVOS REMOVE COMMENT
    write(*,*) "get_lines_tri: lines(1)%normal:",lines(1)%normal,"."
    write(*,*) "get_lines_tri: lines(2)%normal:",lines(2)%normal,"."
    write(*,*) "get_lines_tri: lines(3)%normal:",lines(3)%normal,"."

    ! And calibrate what is the zero of this lines by dotting with
    ! a point we know to be on it
!    do i=1,3
!      lines(i)%c = dot_product(tri%V(:, i), lines(i)%normal)
!    end do

    lines(1)%c = (-1) * lines(1)%normal(1) * tri%V(1,1) - lines(1)%normal(2) * tri%V(2,1)
    lines(2)%c = (-1) * lines(2)%normal(1) * tri%V(1,2) - lines(2)%normal(2) * tri%V(2,2)
    lines(3)%c = (-1) * lines(3)%normal(1) * tri%V(1,3) - lines(3)%normal(2) * tri%V(2,3)
    
    ! IAKOVOS REMOVE COMMENT
    write(*,*) "get_lines_tri: lines(1)%c:",lines(1)%c,"."
    write(*,*) "get_lines_tri: lines(2)%c:",lines(2)%c,"."
    write(*,*) "get_lines_tri: lines(3)%c:",lines(3)%c,"."
    
  end function get_lines_tri

!  pure function unit_cross(vecA, vecB) result(cross)
  pure function unit_cross(vecA) result(cross)
    real, dimension(2), intent(in) :: vecA
    real, dimension(2) :: cross
    cross(1) = vecA(2)
    cross(2) = vecA(1)
!    cross(2) = vecA(3) * vecB(1) - vecA(1) * vecB(3)
!    cross(3) = vecA(1) * vecB(2) - vecA(2) * vecB(1)

!    cross = cross / norm2(cross)
  end function unit_cross

  pure function distances_to_line(line, tri) result(dists)
    type(line_type), intent(in) :: line
    type(tri_type), intent(in) :: tri
    real, dimension(3) :: dists
    integer :: i
    
    forall(i=1:3)
      dists(i) = dot_product(line%normal, tri%V(:, i)) + line%c
    end forall

  end function distances_to_line

  function face_no(i, j, k) result(face)
    ! Given three local node numbers, what is the face that they share?
    integer, intent(in) :: i, j, k
    integer :: face

    do face=1,4
      if (face /= i .and. face /= j .and. face /= k) return
    end do

  end function face_no

end module libsupermesh_tri_intersection_module
