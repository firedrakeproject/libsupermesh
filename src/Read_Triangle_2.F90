module libsupermesh_read_triangle_2

  use libsupermesh_fields_dummy
  use libsupermesh_futils, only : free_unit

  implicit none

  private

  public :: read_triangle_files

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

    unit = free_unit()
    open(unit = unit, file = trim(filename), status = "old", action = "read")

    read(unit, *) nnodes, ldim, nattrs, nbm
    if(nnodes < 0) then
      write(0, "(a)") "Invalid number of vertices"
      stop 1
    end if
    if(ldim /= dim) then
      write(0, "(a)") "Invalid dimension"
      stop 1
    end if
    if(nattrs < 0) then
      write(0, "(a)") "Invalid number of attributes"
      stop 1
    end if
    if(nbm < 0 .or. nbm > 1) then
      write(0, "(a)") "Invalid number of boundary markers"
      stop 1
    end if

    allocate(coords(dim, nnodes))
    if(present(attributes)) allocate(attributes(nattrs, nnodes))
    if(present(boundary_markers)) allocate(boundary_markers(nbm, nnodes))

    allocate(attribute(nattrs), boundary_marker(nbm))
    do i = 1, nnodes
      read(unit, *) ind, coord, attribute, boundary_marker
      if(i /= ind) then
        write(0, "(a)") "Invalid vertex number"
        stop 1
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
  ! Expected format differs from that in the Triangle documentation, in that
  ! arbitrary dimensional linear simplices are supported, and quadratic
  ! simplices are not supported.
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

    unit = free_unit()
    open(unit = unit, file = trim(filename), status = "old", action = "read")

    read(unit, *) nelements, ncell_nodes, nattrs
    if(nelements < 0) then
      write(0, "(a)") "Invalid number of simplices"
      stop 1
    end if
    if(ncell_nodes /= dim + 1) then
      write(0, "(a)") "Invalid number of nodes per simplex"
      stop 1
    end if
    if(nattrs < 0) then
      write(0, "(a)") "Invalid number of attributes"
      stop 1
    end if

    allocate(enlist(ncell_nodes, nelements))
    if(present(attributes)) allocate(attributes(nattrs, nelements))

    allocate(cell_nodes(ncell_nodes), attribute(nattrs))
    do i = 1, nelements
      read(unit, *) ind, cell_nodes, attribute
      if(i /= ind) then
        write(0, "(a)") "Invalid simplex number"
        stop 1
      end if
      if(any(cell_nodes < 1)) then
        write(0, "(a)") "Invalid node"
        stop 1
      end if
      if(present(nnodes)) then
        if(any(cell_nodes > nnodes)) then
          write(0, "(a)") "Invalid node"
          stop 1
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

end module libsupermesh_read_triangle_2
