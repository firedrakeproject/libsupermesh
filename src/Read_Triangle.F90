#include "fdebug.h"

module libsupermesh_read_triangle

  use libsupermesh_debug
  use libsupermesh_fields

  implicit none

  private

  public :: read_node, read_ele, read_triangle_files

contains

  ! Read Triangle .node file, as described at:
  !   https://www.cs.cmu.edu/~quake/triangle.node.html
  ! Comments and empty lines are not supported, and vertices must be indexed
  ! from one. Trailing lines are ignored.
  !
  ! Expected format differs from that in the Triangle documentation, in that
  ! arbitrary dimensional nodal data is supported.
  subroutine read_node(filename, dim, coords, attributes, boundary_markers)
    character(len = *), intent(in) :: filename
    integer, intent(in) :: dim
    ! dim x nnodes
    real, dimension(:, :), allocatable, intent(out) :: coords
    ! nattrs x nnodes
    real, dimension(:, :), allocatable, optional, intent(out) :: attributes
    ! nbm x nnodes
    integer, dimension(:, :), allocatable, optional, intent(out) :: boundary_markers

    integer :: i, ind, nnodes, ldim, nattrs, nbm, unit
    integer, dimension(:), allocatable :: boundary_marker
    real, dimension(dim) :: coord
    real, dimension(:), allocatable :: attribute

    open(newunit = unit, file = trim(filename), status = "old", action = "read")

    read(unit, *) nnodes, ldim, nattrs, nbm
    if(nnodes < 0) then
      libsupermesh_abort("Invalid number of vertices")
    end if
    if(ldim /= dim) then
      libsupermesh_abort("Invalid dimension")
    end if
    if(nattrs < 0) then
      libsupermesh_abort("Invalid number of attributes")
    end if
    if(nbm < 0 .or. nbm > 1) then
      libsupermesh_abort("Invalid number of boundary markers")
    end if

    allocate(coords(dim, nnodes))
    if(present(attributes)) allocate(attributes(nattrs, nnodes))
    if(present(boundary_markers)) allocate(boundary_markers(nbm, nnodes))

    allocate(attribute(nattrs), boundary_marker(nbm))
    do i = 1, nnodes
      if(nattrs > 0) then
        if(nbm > 0) then
          read(unit, *) ind, coord, attribute, boundary_marker
        else
          read(unit, *) ind, coord, attribute
        end if
      else if(nbm > 0) then
        read(unit, *) ind, coord, boundary_marker
      else
        read(unit, *) ind, coord
      end if
      if(i /= ind) then
        libsupermesh_abort("Invalid vertex number")
      end if
      coords(:, i) = coord
      if(present(attributes) .and. nattrs > 0) attributes(:, i) = attribute
      if(present(boundary_markers) .and. nbm > 0) boundary_markers(:, i) = boundary_marker
    end do
    deallocate(attribute, boundary_marker)

    close(unit)

  end subroutine read_node

  ! Read Triangle .ele file, as described at:
  !   https://www.cs.cmu.edu/~quake/triangle.ele.html
  ! Comments and empty lines are not supported, and triangles must be indexed
  ! from one. Trailing lines are ignored.
  !
  ! Expected format differs from that in the Triangle documentation, in that an
  ! arbitrary number of local nodes is permitted.
  subroutine read_ele(filename, dim, enlist, attributes, nnodes)
    character(len = *), intent(in) :: filename
    integer, intent(in) :: dim
    ! loc x nelements
    integer, dimension(:, :), allocatable, intent(out) :: enlist
    ! nattrs x nelements
    real, dimension(:, :), allocatable, optional, intent(out) :: attributes
    integer, optional, intent(in) :: nnodes

    integer :: i, ind, ncell_nodes, nelements, nattrs, unit
    integer, dimension(:), allocatable :: cell_nodes
    real, dimension(:), allocatable :: attribute

    open(newunit = unit, file = trim(filename), status = "old", action = "read")

    read(unit, *) nelements, ncell_nodes, nattrs
    if(nelements < 0) then
      libsupermesh_abort("Invalid number of elements")
    end if
    if(nattrs < 0) then
      libsupermesh_abort("Invalid number of attributes")
    end if

    allocate(enlist(ncell_nodes, nelements))
    if(present(attributes)) allocate(attributes(nattrs, nelements))

    allocate(cell_nodes(ncell_nodes), attribute(nattrs))
    do i = 1, nelements
      if(nattrs > 0) then
        read(unit, *) ind, cell_nodes, attribute
      else
        read(unit, *) ind, cell_nodes
      end if
      if(i /= ind) then
        libsupermesh_abort("Invalid element number")
      end if
      if(any(cell_nodes < 1)) then
        libsupermesh_abort("Invalid node")
      end if
      if(present(nnodes)) then
        if(any(cell_nodes > nnodes)) then
          libsupermesh_abort("Invalid node")
        end if
      end if
      enlist(:, i) = cell_nodes
      if(present(attributes) .and. nattrs > 0) attributes(:, i) = attribute
    end do
    deallocate(cell_nodes, attribute)

    close(unit)

  end subroutine read_ele

  function read_triangle_files(filename, dim) result(positions)
    character(len = *), intent(in) :: filename
    integer, intent(in) :: dim
    
    type(vector_field) :: positions

    integer :: loc, nelements, nnodes
    integer, dimension(:, :), allocatable :: enlist
    real, dimension(:, :), allocatable :: coords
    type(mesh_type) :: mesh

    call read_node(trim(filename) // ".node", dim, coords)
    nnodes = size(coords, 2)
    call read_ele(trim(filename) // ".ele", dim, enlist, nnodes = nnodes)
    loc = size(enlist, 1)
    nelements = size(enlist, 2)

    call allocate(mesh, dim, nnodes, nelements, loc)
    mesh%ndglno = reshape(enlist, (/loc * nelements/))
    deallocate(enlist)
    
    call allocate(positions, dim, mesh)
    positions%val = coords
    deallocate(coords)
    call deallocate(mesh)

  end function read_triangle_files

end module libsupermesh_read_triangle
