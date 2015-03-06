#define BUF_SIZE 150
#include "fdebug.h"

module libsupermesh_tet_intersection_module

  use libsupermesh_elements
!  use vector_tools		! IAKOVOS commented out
  use libsupermesh_fields_data_types
  use libsupermesh_fields_base
  use libsupermesh_fields_allocates
  use libsupermesh_fields_manipulation
  use libsupermesh_transform_elements
  implicit none

  type tet_type
    real, dimension(3, 4) :: V ! vertices of the tet
    integer, dimension(4) :: colours = -1 ! surface colours
  end type tet_type

  type plane_type
    real, dimension(3) :: normal
    real :: c
  end type plane_type

  type(tet_type), dimension(BUF_SIZE), save :: tet_array, tet_array_tmp
  integer :: tet_cnt = 0, tet_cnt_tmp = 0
  type(mesh_type), save :: intersection_mesh
  logical, save :: mesh_allocated = .false.

! IAKOVOS commented out
!  public :: tet_type, plane_type, intersect_tets, get_planes, finalise_tet_intersector
  public :: libsupermesh_intersect_tets

  interface libsupermesh_intersect_tets
    module procedure libsupermesh_intersect_tets_dt
  end interface

  contains

  subroutine finalise_tet_intersector
    if (mesh_allocated) then
!      call deallocate(intersection_mesh)		! ToDo 1
      mesh_allocated = .false.
    end if
  end subroutine finalise_tet_intersector

  subroutine libsupermesh_intersect_tets_dt(tetA, planesB, shape, stat, output, surface_shape, surface_positions, surface_colours)
    type(tet_type), intent(in) :: tetA
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

    tet_cnt = 1
    tet_array(1) = tetA

    if (.not. mesh_allocated) then
      call allocate(intersection_mesh, BUF_SIZE * 4, BUF_SIZE, shape, name="IntersectionMesh")
      intersection_mesh%ndglno = (/ (i, i=1,BUF_SIZE*4) /)
      intersection_mesh%continuity = -1
      mesh_allocated = .true.
    end if

    do i=1,size(planesB)
      ! Clip the tet_array against the i'th plane
      tet_cnt_tmp = 0

      do j=1,tet_cnt
        call clip(planesB(i), tet_array(j))
      end do

      if (i /= size(planesB)) then
        tet_cnt = tet_cnt_tmp
        tet_array(1:tet_cnt) = tet_array_tmp(1:tet_cnt)
      else
        ! Copy the result if the volume is > epsilon
        tet_cnt = 0
        do j=1,tet_cnt_tmp
          vol = tet_volume(tet_array_tmp(j))
          if (vol < 0.0) then
            vec_tmp = tet_array_tmp(j)%V(:, 1)
            colour_tmp = tet_array_tmp(j)%colours(1)
            tet_array_tmp(j)%V(:, 1) = tet_array_tmp(j)%V(:, 2)
            tet_array_tmp(j)%colours(1) = tet_array_tmp(j)%colours(2)
            tet_array_tmp(j)%V(:, 2) = vec_tmp
            tet_array_tmp(j)%colours(2) = colour_tmp
            vol = -vol
          end if

          if (vol > epsilon(0.0)) then
            tet_cnt = tet_cnt + 1
            tet_array(tet_cnt) = tet_array_tmp(j)
          end if
        end do
      end if
    end do

    if (tet_cnt == 0) then
      stat=1
      return
    end if

    stat = 0
    intersection_mesh%nodes = tet_cnt*4
    intersection_mesh%elements = tet_cnt
!    call allocate(output, 3, intersection_mesh, "IntersectionCoordinates")	! ToDo 1

    do ele=1,tet_cnt
      call set(output, ele_nodes(output, ele), tet_array(ele)%V)
    end do

    if (present(surface_positions)) then
      ! OK! Let's loop through all the tets we have and see which faces have positive
      ! colour. These are the ones we want to record in the mesh
      surface_eles = 0
      do ele=1,tet_cnt
        surface_eles = surface_eles + count(tet_array(ele)%colours > 0)
      end do

      call allocate(surface_mesh, surface_eles * 3, surface_eles, surface_shape, name="SurfaceMesh")
      surface_mesh%ndglno = (/ (i, i=1,surface_eles * 3) /)
!      call allocate(surface_positions, 3, surface_mesh, "OutputSurfaceCoordinate")	! ToDo 1
      pwc_surface_mesh = piecewise_constant_mesh(surface_mesh, "PWCSurfaceMesh")
!      call allocate(surface_colours, pwc_surface_mesh, "SurfaceColours")		! ToDo 1
!      call deallocate(surface_mesh)							! ToDo 1
!      call deallocate(pwc_surface_mesh)						! ToDo 1

      j = 1
      do ele=1,tet_cnt
        do i=1,4
          if (tet_array(ele)%colours(i) > 0) then

            ! In python, this is
            ! idx_tmp = [x for x in range(4) if x != i]
            ! Hopefully that will make it clearer
            k = 1
            do l=1,4
              if (l /= i) then
                idx_tmp(k) = l
                k = k + 1
              end if
            end do
            call set(surface_positions, ele_nodes(surface_positions, j), tet_array(ele)%V(:, idx_tmp))
            call set(surface_colours, j, float(tet_array(ele)%colours(i)))
            j = j + 1
          end if
        end do
      end do
    end if

  end subroutine libsupermesh_intersect_tets_dt
  
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

    vol = dot_product(vecA, cross) / 6.0
  end function tet_volume

end module libsupermesh_tet_intersection_module
