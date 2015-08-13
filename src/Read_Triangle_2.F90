#include "fdebug.h"

module libsupermesh_read_triangle_2

  use libsupermesh_fields_dummy
  use libsupermesh_fldebug

  implicit none

  private

  public :: read_triangle_files, read_halo_files

contains
  function free_unit()
    !!< Find a free unit number. Start from unit 10 in order to ensure that
    !!< we skip any preconnected units which may not be correctly identified
    !!< on some compilers.
    integer :: free_unit

    logical :: connected

    do free_unit=10, 99

       inquire(unit=free_unit, opened=connected)

       if (.not.connected) return

    end do

    FLAbort("No free unit numbers avalable")

  end function

  ! Read Triangle .node file, as described at:
  !   https://www.cs.cmu.edu/~quake/triangle.node.html
  ! Comments and empty lines are not supported, and vertices must be indexed
  ! from one. Trailing lines are ignored.
  !
  ! Expected format differs from that in the Triangle documentation, in that
  ! arbitrary dimensional nodal data is supported.
  subroutine read_node(filename, dim, coords, attributes, boundary_markers, nodes)
    character(len = *), intent(in) :: filename
    integer, intent(in) :: dim
    ! dim x nnodes
    real, dimension(:, :), allocatable, intent(out) :: coords
    ! nattrs x nnodes
    real, dimension(:, :), allocatable, optional, intent(out) :: attributes
    ! nbm x nnodes
    integer, dimension(:, :), allocatable, optional, intent(out) :: boundary_markers
    ! number of nodes to read
    integer, optional, intent(in) :: nodes

    integer :: i, ind, nnodes, ldim, nattrs, nbm, unit
    integer, dimension(:), allocatable :: boundary_marker
    real, dimension(dim) :: coord
    real, dimension(:), allocatable :: attribute

    unit = free_unit()
    open(unit = unit, file = trim(filename), status = "old", action = "read")

    read(unit, *) nnodes, ldim, nattrs, nbm
    if(nnodes < 0) then
      FLAbort("Invalid number of vertices")
    end if
    if(ldim /= dim) then
      FLAbort("Invalid dimension")
    end if
    if(nattrs < 0) then
      FLAbort("Invalid number of attributes")
    end if
    if(nbm < 0 .or. nbm > 1) then
      FLAbort("Invalid number of boundary markers")
    end if

    if(present(nodes)) nnodes=nodes
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
        FLAbort("Invalid vertex number")
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

    integer :: i, ind, ncell_nodes, nelements, nattrs, unit, my_elements
    integer, dimension(:), allocatable :: cell_nodes
    real, dimension(:), allocatable :: attribute

    integer, dimension(:, :), allocatable :: enlist_temp
    real, dimension(:, :), allocatable    :: attributes_temp

    unit = free_unit()
    open(unit = unit, file = trim(filename), status = "old", action = "read")

    read(unit, *) nelements, ncell_nodes, nattrs
    if(nelements < 0) then
      FLAbort("Invalid number of simplices")
    end if
    if(ncell_nodes /= dim + 1) then
      FLAbort( "Invalid number of nodes per simplex")
    end if
    if(nattrs < 0) then
      FLAbort("Invalid number of attributes")
    end if

    if(present(nnodes)) then
      allocate(enlist_temp(ncell_nodes,nnodes*10))
    else
      allocate(enlist_temp(ncell_nodes,nelements))
    end if
!    allocate(enlist(ncell_nodes, nelements))
!    if(present(attributes)) allocate(attributes(nattrs, nelements))
    if(present(attributes)) then
      if(present(nnodes)) then
        allocate(attributes_temp(nattrs,nnodes*10))
      else
        allocate(attributes_temp(nattrs,nelements))
      end if
    end if

    enlist_temp = 0
    if(present(attributes)) attributes_temp = 0.0

    allocate(cell_nodes(ncell_nodes), attribute(nattrs))
    my_elements = 0
    do i = 1, nelements
      if(nattrs > 0) then
        read(unit, *) ind, cell_nodes, attribute
      else
        read(unit, *) ind, cell_nodes
      end if
      if(i /= ind) then
        FLAbort("Invalid simplex number")
      end if
      if(any(cell_nodes < 1)) then
        FLAbort("Invalid node")
      end if
      if(present(nnodes)) then
        if(any(cell_nodes > nnodes)) then
!!          FLAbort("Invalid node")
           cycle
        end if
      end if
      my_elements = my_elements + 1
      enlist_temp(:, my_elements) = cell_nodes
!  ToDo TODO todo FIX HACK
      write (*,*) "read_ele: i:",i,", my_elements:",my_elements,", cell_nodes:", cell_nodes
      if(present(attributes) .and. nattrs > 0) attributes_temp(:, my_elements) = attribute
    end do
    deallocate(cell_nodes, attribute)
    
    allocate(enlist(ncell_nodes, my_elements))
    if(present(attributes)) allocate(attributes(nattrs, my_elements))
    do i = 1, my_elements
      enlist(:, i) = enlist_temp(:, i)
      if(present(attributes) .and. nattrs > 0) attributes(:, i) = attributes_temp(:, i)
    end do
    
    deallocate(enlist_temp)
    if(present(attributes)) deallocate(attributes_temp)
!  ToDo TODO todo FIX HACK
    print "(a,I5,a,I5,a)", "read_ele: my_elements:", my_elements,", all elements:",nelements,"."

    close(unit)

  end subroutine read_ele
  
  subroutine read_halo(filename, nodes)
    character(len = *), intent(in) :: filename
    integer, intent(out) :: nodes

    character(len=1000) :: buffer
    integer :: unit, ios = 0, line = 0, pos = 0
  
    unit = free_unit()
    open(unit = unit, file = trim(filename), status = "old", action = "read")

    do while (ios == 0)
      read(unit, '(A)', iostat=ios) buffer
      if (ios == 0) then
        line = line + 1
        pos = 0

        ! Find the first instance of n_private_nodes=.
        pos = index(buffer, 'n_private_nodes')

        if ( pos .gt. 1 ) then
          buffer = buffer(pos+1+16:)
          pos = index(buffer, '"')
          buffer = buffer(:pos-1)
          read (buffer,'(I10)') nodes
        end if

      end if
    end do

  end subroutine read_halo

  function read_triangle_files(filename, dim, nodes) result(positions)
    character(len = *), intent(in) :: filename
    integer, intent(in) :: dim
    integer, optional, intent(in) :: nodes
    
    type(vector_field) :: positions

    integer :: loc, nelements, nnodes
    integer, dimension(:, :), allocatable :: enlist
    real, dimension(:, :), allocatable :: coords
    type(mesh_type) :: mesh

    call read_node(trim(filename) // ".node", dim, coords, nodes = nodes)
    nnodes = size(coords, 2)
!  ToDo TODO todo FIX HACK
    print "(a,I5)", "read_triangle_files: nnodes (cells):", nnodes
    call read_ele(trim(filename) // ".ele", dim, enlist, nnodes = nnodes)
    loc = size(enlist, 1)
    nelements = size(enlist, 2)
!  ToDo TODO todo FIX HACK
    print "(a,I5)", "read_triangle_files: elements (tris):", nelements

    call allocate(mesh, dim, nnodes, nelements, loc)
    mesh%ndglno = reshape(enlist, (/loc * nelements/))
    deallocate(enlist)
    
    call allocate(positions, dim, mesh)
    positions%val = coords
    deallocate(coords)
    call deallocate(mesh)
    
  end function read_triangle_files
  
  function read_halo_files(filename) result(nodes)
    character(len = *), intent(in) :: filename
    integer :: nodes
    
    call read_halo(trim(filename) // ".halo", nodes)
  end function read_halo_files

end module libsupermesh_read_triangle_2
