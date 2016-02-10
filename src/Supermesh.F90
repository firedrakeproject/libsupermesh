#include "libsupermesh_debug.h"

module libsupermesh_supermesh

  use libsupermesh_debug, only : abort_pinpoint
  use libsupermesh_interval_intersection, only : intersect_intervals, &
    & interval_buf_size
  use libsupermesh_tri_intersection, only : tri_type, line_type, max_n_trisC, &
    & intersect_tris, tri_buf_size, intersect_polys, get_lines, triangle_area
  use libsupermesh_tet_intersection, only : tet_type, plane_type, max_n_tetsC, &
    & intersect_tets, tet_buf_size, tetrahedron_volume

  implicit none

  private

  public :: max_n_simplicesC, intersection_buffer_size, intersect_simplices, &
    & intersect_elements, simplex_volume, triangle_area, tetrahedron_volume

contains

  pure function max_n_simplicesC(dim) result(size)
    integer, intent(in) :: dim

    integer :: size

    select case(dim)
      case(1)
        size = interval_buf_size
      case(2)
        size = tri_buf_size
      case(3)
        size = tet_buf_size
      case default
        size = 0
    end select

  end function max_n_simplicesC

  pure function intersection_buffer_size(dim, loc_a, loc_b) result(size)
    integer, intent(in) :: dim
    integer, intent(in) :: loc_a
    integer, intent(in) :: loc_b

    integer :: size

    if(loc_a == dim + 1 .and. loc_b == dim + 1) then
      size = max_n_simplicesC(dim)
    else if((loc_a == dim + 1 .and. loc_b == 2 ** dim) &
     & .or. (loc_b == dim + 1 .and. loc_a == 2 ** dim)) then
      select case(dim)
        case(2)
          size = max_n_trisC(n_linesB = 4)
        case(3)
          size = max_n_tetsC(n_planesB = 6)
        case default
          size = 0
      end select
    else if(loc_a == 2 ** dim .and. loc_b == 2 ** dim) then
      select case(dim)
        case(2)
          size = max_n_trisC(n_linesA = 4, n_linesB = 4)
        case(3)
          size = 5 * max_n_tetsC(n_planesB = 6)
        case default
          size = 0
      end select
    else
      size = 0
    end if

  end function intersection_buffer_size

  subroutine intersect_simplices(simplexA, simplexB, simplicesC, n_simplicesC)
    ! dim x loc_a
    real, dimension(:, :), intent(in) :: simplexA
    ! dim x loc_b
    real, dimension(:, :), intent(in) :: simplexB
    real, dimension(:, :, :), intent(inout) :: simplicesC
    integer, intent(out) :: n_simplicesC
    
    select case(size(simplexA, 1))
      case(1)
        call intersect_intervals(simplexA(1, :), simplexB(1, :), simplicesC(1, :, 1), n_simplicesC)
      case(2)
        call intersect_tris(simplexA, simplexB, simplicesC, n_simplicesC)
      case(3)
        call intersect_tets(simplexA, simplexB, simplicesC, n_simplicesC)
      case default
        libsupermesh_abort("Unsupported element type")
    end select
    
  end subroutine intersect_simplices

  subroutine intersect_elements(elementA, elementB, elementsC, n_elementsC)
    ! dim x loc_a
    real, dimension(:, :), intent(in) :: elementA
    ! dim x loc_b
    real, dimension(:, :), intent(in) :: elementB
    real, dimension(:, :, :), intent(inout) :: elementsC
    integer, intent(out) :: n_elementsC

    integer :: dim, loc_a, loc_b
  
    dim = size(elementA, 1)
    loc_a = size(elementA, 2)
    loc_b = size(elementB, 2)
    
    if(loc_a == dim + 1 .and. loc_b == dim + 1) then
      call intersect_simplices(elementA, elementB, elementsC, n_elementsC)
    else if(loc_a == dim + 1 .and. loc_b == 2 ** dim) then
      select case(dim)
        case(2)
          call intersect_tri_quad(elementA, elementB, elementsC, n_elementsC)
        case(3)
          call intersect_tet_hex(elementA, elementB, elementsC, n_elementsC)
        case default
          libsupermesh_abort("Unsupported element type")
      end select
    else if(loc_a == 2 ** dim .and. loc_b == dim + 1) then
      select case(dim)
        case(2)
          call intersect_tri_quad(elementB, elementA, elementsC, n_elementsC)
        case(3)
          call intersect_tet_hex(elementB, elementA, elementsC, n_elementsC)
        case default
          libsupermesh_abort("Unsupported element type")
      end select
    else if(loc_a == 2 ** dim .and. loc_b == 2 ** dim) then
      select case(dim)
        case(2)
          call intersect_quads(elementA, elementB, elementsC, n_elementsC)
        case(3)
          call intersect_hexes(elementA, elementB, elementsC, n_elementsC)
        case default
          libsupermesh_abort("Unsupported element type")
      end select
    else
      libsupermesh_abort("Unsupported element type")
    end if

  end subroutine intersect_elements
  
  subroutine intersect_tri_quad(triA, quadB, trisC, n_trisC)
    ! 2 x 3
    real, dimension(:, :), intent(in) :: triA
    ! 2 x 4
    real, dimension(:, :), intent(in) :: quadB
    real, dimension(:, :, :), intent(inout) :: trisC
    integer, intent(out) :: n_trisC
    
    integer :: i
    real, dimension(2, 3 * (2 ** 4), 2), save :: work
    type(tri_type) :: ltriA
    type(tri_type), dimension(3 * (2 ** 4) - 2), save :: ltrisC
    
    ltriA%v = triA
    call intersect_polys(ltriA, get_lines(quadB), ltrisC, n_trisC, work = work)
    do i = 1, n_trisC
      trisC(:, :, i) = ltrisC(i)%v
    end do
  
  end subroutine intersect_tri_quad

  subroutine intersect_quads(quadA, quadB, trisC, n_trisC)
    ! 2 x 4
    real, dimension(:, :), intent(in) :: quadA
    ! 2 x 4
    real, dimension(:, :), intent(in) :: quadB
    real, dimension(:, :, :), intent(inout) :: trisC
    integer, intent(out) :: n_trisC
    
    integer :: i
    real, dimension(2, 4 * (2 ** 4), 2), save :: work
    type(tri_type), dimension(4 * (2 ** 4) - 2), save :: ltrisC
    
    call intersect_polys(quadA, quadB, ltrisC, n_trisC, work = work)
    do i = 1, n_trisC
      trisC(:, :, i) = ltrisC(i)%v
    end do
  
  end subroutine intersect_quads

  subroutine intersect_tet_hex(tetA, hexB, tetsC, n_tetsC)
    ! 3 x 4
    real, dimension(:, :), intent(in) :: tetA
    ! 4 x 8
    real, dimension(:, :), intent(in) :: hexB
    real, dimension(:, :, :), intent(inout) :: tetsC
    integer, intent(out) :: n_tetsC
    
    integer :: i, ln_tetsC
    type(plane_type), dimension(6) :: planesB
    type(tet_type) :: ltetA
    type(tet_type), dimension(3 ** 6), save :: ltetsC, work
    
    ! Gmsh node ordering. See:
    !   http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
    
    ltetA%v = tetA
    planesB(1)%normal = unit_cross(hexB(:, 2) - hexB(:, 1), hexB(:, 6) - hexB(:, 1))
    planesB(1)%c = dot_product(hexB(:, 1), planesB(1)%normal)
    planesB(2)%normal = unit_cross(hexB(:, 6) - hexB(:, 5), hexB(:, 7) - hexB(:, 5))
    planesB(2)%c = dot_product(hexB(:, 5), planesB(2)%normal)
    planesB(3)%normal = unit_cross(hexB(:, 2) - hexB(:, 6), hexB(:, 3) - hexB(:, 6))
    planesB(3)%c = dot_product(hexB(:, 6), planesB(3)%normal)
    planesB(4)%normal = unit_cross(hexB(:, 1) - hexB(:, 2), hexB(:, 4) - hexB(:, 2))
    planesB(4)%c = dot_product(hexB(:, 2), planesB(4)%normal)
    planesB(5)%normal = unit_cross(hexB(:, 5) - hexB(:, 1), hexB(:, 8) - hexB(:, 1))
    planesB(5)%c = dot_product(hexB(:, 1), planesB(5)%normal)
    planesB(6)%normal = unit_cross(hexB(:, 7) - hexB(:, 8), hexB(:, 3) - hexB(:, 8))
    planesB(6)%c = dot_product(hexB(:, 8), planesB(6)%normal)
    
    call intersect_tets(ltetA, planesB, ltetsC, ln_tetsC, work = work)
    do i = 1, ln_tetsC
      tetsC(:, :, i) = ltetsC(i)%v
    end do
    n_tetsC = ln_tetsC
   
  contains

    pure function unit_cross(vecA, vecB) result(cross)
      real, dimension(3), intent(in) :: vecA, vecB
      real, dimension(3) :: cross
      
      cross(1) = vecA(2) * vecB(3) - vecA(3) * vecB(2)
      cross(2) = vecA(3) * vecB(1) - vecA(1) * vecB(3)
      cross(3) = vecA(1) * vecB(2) - vecA(2) * vecB(1)
      cross = cross / norm2(cross)
      
    end function unit_cross
   
  end subroutine intersect_tet_hex

  subroutine intersect_hexes(hexA, hexB, tetsC, n_tetsC)
    ! 3 x 8
    real, dimension(:, :), intent(in) :: hexA
    ! 3 x 8
    real, dimension(:, :), intent(in) :: hexB
    real, dimension(:, :, :), intent(inout) :: tetsC
    integer, intent(out) :: n_tetsC
    
    integer :: ln_tetsC
    
    ! Gmsh node ordering. See:
    !   http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
    
    ! Cube slicing as, e.g. in:
    !   Isosurfaces: Geometry, Topology, and Algorithms, R. Wenger, CRC Press,
    !   2013, figure 2.29
    
    ! Slicing off two opposite corners on the top ...
    call intersect_tet_hex(hexA(:, (/3, 6, 7, 8/)), hexB, tetsC, ln_tetsC)
    n_tetsC = ln_tetsC
    call intersect_tet_hex(hexA(:, (/1, 3, 4, 8/)), hexB, tetsC(:, :, 1 + n_tetsC:), ln_tetsC)
    n_tetsC = n_tetsC + ln_tetsC
    ! ... and two opposite corners on the bottom ...
    call intersect_tet_hex(hexA(:, (/1, 5, 6, 8/)), hexB, tetsC(:, :, 1 + n_tetsC:), ln_tetsC)
    n_tetsC = n_tetsC + ln_tetsC
    call intersect_tet_hex(hexA(:, (/1, 2, 3, 6/)), hexB, tetsC(:, :, 1 + n_tetsC:), ln_tetsC)
    n_tetsC = n_tetsC + ln_tetsC
    ! ... to leave a single tetrahedron in the centre
    call intersect_tet_hex(hexA(:, (/1, 3, 6, 8/)), hexB, tetsC(:, :, 1 + n_tetsC:), ln_tetsC)
    n_tetsC = n_tetsC + ln_tetsC
  
  end subroutine intersect_hexes

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

end module libsupermesh_supermesh
