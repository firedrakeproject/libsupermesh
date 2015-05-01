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
    real :: c = 0
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
        call clip(linesB(i), tri_array(j))
      end do
      
      ! IAKOVOS REMOVE COMMENT
      write(*,*) "libsupermesh_intersect_tris_dt: 3: i:",i,", j:",j,", tri_cnt_tmp:",tri_cnt_tmp,", tri_cnt:",tri_cnt,"."
      
      if (i /= size(linesB)) then
        ! IAKOVOS REMOVE COMMENT
        write(*,*) "libsupermesh_intersect_tris_dt: Inside IF i:",i,"."
        tri_cnt = tri_cnt_tmp
        tri_array(1:tri_cnt) = tri_array_tmp(1:tri_cnt)
      else
        ! IAKOVOS REMOVE COMMENT
        write(*,*) "libsupermesh_intersect_tris_dt: Inside ELSE i:",i,", j:",j,"."
        tri_cnt = tri_cnt
!        tri_array(tri_cnt) = tri_array_tmp(j)
        tri_array(1:tri_cnt) = tri_array_tmp(1:tri_cnt)
      end if
    end do
    
    ! IAKOVOS REMOVE COMMENT
    write(*,*) "libsupermesh_intersect_tris_dt: 4: tri_cnt:",tri_cnt,"."
    
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
      ! IAKOVOS REMOVE COMMENT
      write(*,*) "libsupermesh_intersect_tris_dt: 6: ele_nodes(output, ele):",ele_nodes(output, ele),", tri_array(",ele,")%V:",tri_array(ele)%V,"."
      
      call set(output, ele_nodes(output, ele), tri_array(ele)%V)
    end do
    
    ! IAKOVOS REMOVE COMMENT
    write(*,*) "libsupermesh_intersect_tris_dt: 7"
    
    if (present(surface_positions)) then
      FLAbort("libsupermesh_intersect_tris_dt: Not surface_positions")
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
          w0 = -dists(neg_idx(1)) * invdiff
          w1 =  dists(pos_idx(i)) * invdiff
          
          write(*,*) "clip: Case 1: subcase 2:",i," tri_cnt_tmp:",tri_cnt_tmp,", point:",pos_idx(i),"(",tri%V(:,pos_idx(i)),") is POSITIVE (invdiff:",invdiff,")."
          write(*,*) "clip: Case 1: subcase 2:",i," w0:",w0,",w1:",w1,"."
          write(*,*) "clip: Case 1: subcase 2:",i," X and Y:",w0 * tri_array_tmp(tri_cnt_tmp)%V(:, pos_idx(i)) + w1 * tri_array_tmp(tri_cnt_tmp)%V(:, neg_idx(1)),"."
          
!          tri_array_tmp(tri_cnt_tmp)%V(:, pos_idx(i)) = &
!           w0 * tri_array_tmp(tri_cnt_tmp)%V(:, pos_idx(i)) + &
!           w1 * tri_array_tmp(tri_cnt_tmp)%V(:, neg_idx(1))
        
        end do
        do i=1,neg_cnt
          invdiff = 1.0 / ( dists(pos_idx(1)) - dists(neg_idx(i)) )
          write(*,*) "clip: Case 1: subcase 2:",i," point:",neg_idx(i),"(",tri%V(:,neg_idx(i)),") is NEGATIVE (invdiff:",invdiff,")."
          
          w0 = -dists(neg_idx(i)) * invdiff
          w1 =  dists(pos_idx(1)) * invdiff
          
          write(*,*) "clip: Case 1: subcase 2:",i," w0:",w0,",w1:",w1,"."
          write(*,*) "clip: Case 1: subcase 2:",i," X and Y:",w0 * tri%V(:, pos_idx(1)) + w1 * tri%V(:, neg_idx(i)),"."
          
          tri_array_tmp(tri_cnt_tmp)%V(:, neg_idx(i)) = &
           w0 * tri%V(:, pos_idx(1)) + w1 * tri%V(:, neg_idx(i))
          
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
        
        
        tri_cnt_tmp = tri_cnt_tmp + 1
        tri_array_tmp(tri_cnt_tmp) = tri
        do i=1,pos_cnt
          invdiff = 1.0 / ( dists(pos_idx(i)) - dists(neg_idx(1)) )
          write(*,*) "clip: Case 2: subcase 1:",i," point:",pos_idx(i),"(",tri%V(:,pos_idx(i)),") is POSITIVE (invdiff:",invdiff,")."
!          if ( )
        
        end do
        do i=1,neg_cnt
          invdiff = 1.0 / ( dists(pos_idx(1)) - dists(neg_idx(i)) )
          write(*,*) "clip: Case 2: subcase 1:",i," point:",neg_idx(i),"(",tri%V(:,neg_idx(i)),") is NEGATIVE (invdiff:",invdiff,")."
          
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

    real, dimension(2) :: edge10, edge20, edge30
    real :: det
    integer :: i
    ! P = 1, Q = 2, R = 3
    edge10 = tri%V(:, 2) - tri%V(:, 1); ! PQ = ( Qx - Px , Qy - Py )
    edge20 = tri%V(:, 3) - tri%V(:, 2); ! QR = ( Rx - Qx , Ry - Qy )
    edge30 = tri%V(:, 1) - tri%V(:, 3); ! RP = ( Px - Rx , Py - Ry )
    
    ! IAKOVOS REMOVE COMMENT
    write(*,*) "get_lines_tri: input: point1 (P):",tri%V(:,1),", point2 (Q):",tri%V(:,2),&
            ", point3 (R):",tri%V(:,3),"."
    write(*,*) "get_lines_tri: edge10 (PQ):",edge10,"."
    write(*,*) "get_lines_tri: edge20 (QR):",edge20,"."
    write(*,*) "get_lines_tri: edge30 (PR):",edge30,"."
 
    lines(1)%normal = unit_cross(edge10) ! PQ normal
    lines(2)%normal = unit_cross(edge20) ! QR normal
    lines(3)%normal = unit_cross(edge30) ! PR normal
    
    det = dot_product(edge30, lines(3)%normal)
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
    
    ! IAKOVOS REMOVE COMMENT
    do i=1,3
      write(*,*) "get_lines_tri: lines(",i,")%normal:",lines(i)%normal,", lines%c:",lines(i)%c,"."
    end do
    
  end function get_lines_tri
  
  pure function unit_cross(vecA) result(cross)
    real, dimension(2), intent(in) :: vecA
    real, dimension(2) :: cross
    cross(1) = (-1) * vecA(2)
    cross(2) = vecA(1)

  end function unit_cross

!  pure function distances_to_line(line, tri) result(dists)
  function distances_to_line(line, tri) result(dists)
    type(line_type), intent(in) :: line
    type(tri_type), intent(in) :: tri
    real, dimension(3) :: dists
    real, dimension(2) :: P
    integer :: i
    
!    P(1) = 4
!    P(2) = 4
!    P(1,2) = 4
!    P(2,2) = 4
!    P(1,3) = 4
!    P(2,3) = 4
    
!    forall(i=1:3)
    do i=1,3
!      dists(i) = dot_product(line%normal, tri%V(:, i)) - line%c
!      dists(i) = dot_product((tri%V(:,i) - P), line%normal )
      dists(i) = dot_product(line%normal, tri%V(:, i)) - line%c
      
      ! IAKOVOS REMOVE COMMENT
!      write(*,*) "distances_to_line: i:",i," line%normal:",line%normal,"c:",line%c,", tri%V(:,",i,"):",tri%V(:,i)," - P(",i,"):",P," = ",dists(i),"."
      write(*,*) "distances_to_line: i:",i," line%normal:",line%normal,"c:",line%c,", tri%V(:,",i,"):",tri%V(:,i)," = ",dists(i)
!    end forall
    end do

  end function distances_to_line

end module libsupermesh_tri_intersection_module
