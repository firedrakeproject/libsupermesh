#include "fdebug.h"

module libsupermesh_read_triangle

  use libsupermesh_fields
  use libsupermesh_fldebug

  implicit none

  private

  public :: read_triangle_files, read_halo_files

contains

  ! Read Triangle .node file, as described at:
  !   https://www.cs.cmu.edu/~quake/triangle.node.html
  ! Comments and empty lines are not supported, and vertices must be indexed
  ! from one. Trailing lines are ignored.
  !
  ! Expected format differs from that in the Triangle documentation, in that
  ! arbitrary dimensional nodal data is supported.
  subroutine read_node(filename, dim, coords, attributes, boundary_markers, nnodes, nodes, nodes_translate)
    character(len = *), intent(in) :: filename
    integer, intent(in) :: dim
    ! dim x nnodes
    real, dimension(:, :), allocatable, intent(out) :: coords
    ! nattrs x nnodes
    real, dimension(:, :), allocatable, optional, intent(out) :: attributes
    ! nbm x nnodes
    integer, dimension(:, :), allocatable, optional, intent(out) :: boundary_markers
    ! number of nodes to read
    integer, optional, intent(in) :: nnodes
    integer, dimension(:), optional, intent(in)  :: nodes
    integer, dimension(:), optional, intent(out) :: nodes_translate

    integer :: i, j, ind, nnodes_local, ldim, nattrs, nbm, unit, my_node = 0
    integer, dimension(:), allocatable :: boundary_marker
    real, dimension(dim) :: coord
    real, dimension(:), allocatable :: attribute
    logical                         :: found_match

    open(newunit = unit, file = trim(filename), status = "old", action = "read")

    read(unit, *) nnodes_local, ldim, nattrs, nbm
    if(nnodes_local < 0) then
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

    if(present(nnodes)) nnodes_local=nnodes
    if(present(nodes)) then
      allocate(coords(dim,size(nodes)))
      if(present(attributes)) allocate(attributes(nattrs, size(nodes)))
      if(present(boundary_markers)) allocate(boundary_markers(nbm, size(nodes)))
    else
      allocate(coords(dim, nnodes_local))
      if(present(attributes)) allocate(attributes(nattrs, nnodes_local))
      if(present(boundary_markers)) allocate(boundary_markers(nbm, nnodes_local))
    end if

    allocate(attribute(nattrs), boundary_marker(nbm))
    my_node = 0
    do i = 1, nnodes_local
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
      if(present(nodes)) then
        found_match = .false.
        do j = 1, size(nodes)
          if(i .eq. nodes(j)) then
            found_match = .true.
          end if
        end do

        if (found_match .eqv. .false.) then
!write (*,*) "read_node: CYCLE!!! i:",i,", my_node:",my_node,", coord:", coord
          cycle
        end if
      end if

      my_node = my_node + 1
!write (*,*) "read_node: i:",i,", my_node:",my_node,", coord:", coord
      if(present(nodes)) then
        coords(:, my_node) = coord
        nodes_translate(my_node) = i
      else
        coords(:, i) = coord
      end if

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
  subroutine read_ele(filename, dim, enlist, attributes, nnodes, nodes)
    character(len = *), intent(in) :: filename
    integer, intent(in) :: dim
    ! loc x nelements
    integer, dimension(:, :), allocatable, intent(out) :: enlist
    ! nattrs x nelements
    real, dimension(:, :), allocatable, optional, intent(out) :: attributes
    integer, optional, intent(in) :: nnodes
    integer, dimension(:), optional, intent(in) :: nodes

    integer :: i, j, k, ind, ncell_nodes, nelements, nattrs, unit, my_elements
    integer, dimension(:), allocatable :: cell_nodes
    real, dimension(:), allocatable :: attribute

    integer, dimension(:, :), allocatable :: enlist_temp
    real, dimension(:, :), allocatable    :: attributes_temp
    logical, dimension(:), allocatable    :: found_match

    open(newunit = unit, file = trim(filename), status = "old", action = "read")

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

    if(present(nodes)) then
      allocate(enlist_temp(ncell_nodes,size(nodes)*100))
    else
      allocate(enlist_temp(ncell_nodes,nelements))
    end if
!    allocate(enlist(ncell_nodes, nelements))
!    if(present(attributes)) allocate(attributes(nattrs, nelements))
    if(present(attributes)) then
      if(present(nodes)) then
        allocate(attributes_temp(nattrs,size(nodes)*100))
      else
        allocate(attributes_temp(nattrs,nelements))
      end if
    end if

    enlist_temp = 0
    if(present(attributes)) attributes_temp = 0.0

    allocate(cell_nodes(ncell_nodes), attribute(nattrs), found_match(ncell_nodes))
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
!      if(present(nnodes)) then
!        if(any(cell_nodes > nnodes)) then
!          FLAbort("Invalid node.")
!        end if
!      end if
      if(present(nodes)) then
        found_match = .false.
        do j = 1, ncell_nodes
          do k = 1, size(nodes)
            if(cell_nodes(j) .eq. nodes(k)) then
              found_match(j) = .true.
            end if
          end do
        end do

        if (any(found_match .eqv. .false.)) then
!write (*,*) "read_ele: CYCLE!!! i:",i,", my_elements:",my_elements,", cell_nodes:", cell_nodes
call FLUSH()
          cycle
        end if
      end if
      my_elements = my_elements + 1
      enlist_temp(:, my_elements) = cell_nodes
!  ToDo TODO todo FIX HACK
!      write (*,*) "read_ele: i:",i,", my_elements:",my_elements,", cell_nodes:", cell_nodes
      if(present(attributes) .and. nattrs > 0) attributes_temp(:, my_elements) = attribute
    end do
    deallocate(cell_nodes, attribute, found_match)
    
    allocate(enlist(ncell_nodes, my_elements))
    if(present(attributes)) allocate(attributes(nattrs, my_elements))
    do i = 1, my_elements
      enlist(:, i) = enlist_temp(:, i)
      if(present(attributes) .and. nattrs > 0) attributes(:, i) = attributes_temp(:, i)
    end do
    
    deallocate(enlist_temp)
    if(present(attributes)) deallocate(attributes_temp)
!  ToDo TODO todo FIX HACK
!    print "(a,I5,a,I5,a)", "read_ele: my_elements:", my_elements,", all elements:",nelements,"."

    close(unit)

  end subroutine read_ele
  
  subroutine read_halo(filename, nnodes, nodes)
    character(len = *), intent(in) :: filename
    integer, intent(out) :: nnodes
    integer, dimension(:), allocatable, intent(out) :: nodes

    character(len=5000) :: buffer
    integer :: unit, ios, line, pos, halo_data_process_pos, halo_level_end_pos, &
        & halo_data_process_end_pos, halo_data_receive_pos, count_i, count_pos, &
        & count_nwords, temp_value, my_node, other_node, i
    integer, dimension(:), allocatable :: tmp_array
  
    ios = 0; line = 0; pos = 0; halo_data_process_pos = 0; halo_level_end_pos = -1;
    halo_data_process_end_pos = -1; halo_data_receive_pos = 0; count_i = 0; count_pos = 1;
    count_nwords = 0; temp_value = 0; my_node = 0; other_node = 0; i = 0
    open(newunit = unit, file = trim(filename), status = "old", action = "read")

    do while (ios == 0)
      read(unit, '(A)', iostat=ios) buffer
      if (ios == 0) then
        line = line + 1
        pos = 0

        ! Find the first instance of "halos process"
        pos = index(buffer, 'halos process')
        if ( pos .gt. 1 ) then
          buffer = buffer(pos+1+14:)
          pos = index(buffer, '"')
          buffer = buffer(:pos-1)
          read (buffer,'(I10)') my_node
!          print "(a,I5)", "read_halo: my_node:", my_node
        end if

        ! Find the first instance of halo level="1" n_private_nodes=
        pos = index(buffer, '<halo level="1" n_private_nodes=')
        if ( pos .gt. 1 ) then
          buffer = buffer(pos+1+32:)
          pos = index(buffer, '"')
          buffer = buffer(:pos-1)
          read (buffer,'(I10)') nnodes
!          print "(a,I5)", "read_halo: nnodes:", nnodes
!        end if

        ! Find the first instance of halo level="2"
!        pos = index(buffer, 'halo level="2"')
!        if ( pos .gt. 1 ) then
          ! We are now parsing "halo level"s
          DO WHILE ( halo_level_end_pos .lt. 1 )
            read(unit, '(A)', iostat=ios) buffer
            if (ios == 0) then
              line = line + 1
              pos = 0

              halo_level_end_pos = index(buffer, '</halo>')
              if ( halo_level_end_pos .gt. 1 ) cycle

!              print "(a,a)", "read_halo1: buffer:", trim(adjustl(buffer))

              pos = index(buffer, '<halo_data process="')
              if ( pos .gt. 1 ) then
                halo_data_process_end_pos = -1
                DO WHILE ( halo_data_process_end_pos .lt. 1 )
                  ! Check if this is the first line of the tag
                  halo_data_process_pos = index(buffer, '<halo_data process="')
                  if ( halo_data_process_pos .gt. 1 ) then
                    ! If this is the first line of the tag, extract the other node
!                    print "(a,I5,a,a)", "read_halo__: buffer:", halo_data_process_pos, " ", trim(adjustl(buffer))
                    buffer = buffer(halo_data_process_pos+1+19:)
!                    print "(a,a)", "read_halo__: buffer:", trim(adjustl(buffer))
                    halo_data_process_pos = index(buffer, '"')
                    buffer = buffer(:halo_data_process_pos-1)

                    read (buffer,'(I10)') other_node
!                    print "(a,I5)", "read_halo: other_node:", other_node
                  end if

                  halo_data_process_end_pos = index(buffer, '</halo_data>')
                  if ( halo_data_process_end_pos .gt. 1 ) exit  ! Cycle if we found the end of the tag
                  if ( my_node .gt. other_node ) exit           ! Cycle if our node is smaller than the other node

                  read(unit, '(A)', iostat=ios) buffer
                  if (ios == 0) then
                    line = line + 1
                    halo_data_receive_pos = 0

                    ! Find the starting position of the <receive> tag.
                    halo_data_receive_pos = index(buffer, '<receive>')

                    ! If not found cycle?
                    if ( halo_data_receive_pos .lt. 1 ) cycle

                    ! Find the position of the </receive> tag and discard the tags.
                    buffer = buffer(halo_data_receive_pos+9:)
                    halo_data_receive_pos = index(buffer, '</receive>')
                    buffer = buffer(:halo_data_receive_pos-1)

                    ! Count how many nodes we expect to find (in order to allocate/resize the nodes array).
                    count_nwords = count_nodes(buffer)
!                    print "(a,I5)", "read_halo2: count:", count_nwords

                    ! Allocate / Resize the nodes array.
                    if (allocated(nodes) .eqv. .false.) then
                      allocate(nodes(count_nwords + nnodes))
!                      print "(a,I5,a,I5,a)", "read_halo2: count_nwords:", count_nwords,", nnodes:",nnodes,"."
                      count_nwords = nnodes
                      do i = 1, nnodes
                        nodes(i) = i
                      end do
                    else
                      allocate(tmp_array(size(nodes) + count_nwords))
                      tmp_array(1:size(nodes)) = nodes
                      count_nwords = size(nodes)
                      deallocate(nodes)
                      allocate(nodes(size(tmp_array)))
                      nodes = tmp_array
                      deallocate(tmp_array)
!                      print "(a,I5,a,I5,a)", "read_halo2: count_nwords:", count_nwords,", size(nodes):",size(nodes),"."
                    end if

                    ! Get the value of each node
                    call get_nodes(buffer, count_nwords, nodes)
                  end if
                end do
              end if
            end if
          end do
        end if
      end if
    end do

    nnodes = size(nodes)

    close(unit)

  end subroutine read_halo

  function count_nodes(buffer) result(count_nwords)
    character(len = *), intent(in) :: buffer
    integer :: count_nwords

    integer :: count_i, count_pos

    count_i = 0; count_pos = 1; count_nwords = 0
    getnodescount: do
      count_i = verify(buffer(count_pos:), ' ')!-- Find next non-blank.
      if (count_i == 0) exit getnodescount     !-- No word found.
      count_nwords = count_nwords + 1          !-- Found something.
      count_pos = count_pos + count_i - 1      !-- Move to start of the word.

      count_i = scan(buffer(count_pos:), ' ')  !-- Find next blank.
      if (count_i == 0) exit getnodescount     !-- No blank found.
      count_pos = count_pos + count_i - 1      !-- Move to the blank.
    end do getnodescount
  end function count_nodes

  subroutine get_nodes(buffer, count_nwords, nodes)
    character(len = *), intent(in)       :: buffer
    integer, intent(inout)               :: count_nwords
    integer, dimension(:), intent(inout) :: nodes

    integer :: count_i, count_pos, temp_value

    count_i = 0; count_pos = 1
    getnodes: do
      count_i = verify(buffer(count_pos:), ' ')!-- Find next non-blank.
      if (count_i == 0) exit getnodes          !-- No word found.
      count_nwords = count_nwords + 1          !-- Found something.
      count_pos = count_pos + count_i - 1      !-- Move to start of the word.

      count_i = scan(buffer(count_pos:), ' ')  !-- Find next blank.
      if (count_i == 0) exit getnodes          !-- No blank found.
!      print "(a,a,a)", "read_halo2: buffer__:", buffer(count_pos:count_pos + count_i - 2),"."
      read (buffer(count_pos:count_pos + count_i - 2),'(I10)') temp_value
      nodes(count_nwords) = temp_value
      count_pos = count_pos + count_i - 1      !-- Move to the blank.
    end do getnodes
  end subroutine get_nodes

  function read_triangle_files(filename, dim, nnodes, nodes) result(positions)
    character(len = *), intent(in) :: filename
    integer, intent(in) :: dim
    integer, optional, intent(in) :: nnodes
    integer, dimension(:), optional, intent(in) :: nodes
    
    type(vector_field) :: positions

    integer :: loc, nelements, nnodes_local
    integer, dimension(:, :), allocatable :: enlist
    integer, dimension(:), allocatable :: nodes_translate
    real, dimension(:, :), allocatable :: coords
!    integer, dimension(dim+1)  :: cell_nodes
    type(mesh_type) :: mesh

    if(present(nodes)) then
      allocate(nodes_translate(size(nodes)))
      call read_node(trim(filename) // ".node", dim, coords, nodes = nodes, nodes_translate = nodes_translate)
    else
      call read_node(trim(filename) // ".node", dim, coords, nnodes = nnodes)
    end if
    nnodes_local = size(coords, 2)

!  ToDo TODO todo FIX HACK
!    print "(a,I5)", "read_triangle_files: nnodes_local (cells):", nnodes_local
    if (present(nodes)) then
      call read_ele(trim(filename) // ".ele", dim, enlist, nnodes = nnodes_local, nodes = nodes)
    else
      call read_ele(trim(filename) // ".ele", dim, enlist, nnodes = nnodes_local)
    end if
    loc = size(enlist, 1)
    nelements = size(enlist, 2)
!  ToDo TODO todo FIX HACK
!    print "(a,I5)", "read_triangle_files: elements (tris):", nelements

    ! Fix enlist, replacing node/element order and IDs.
    if(present(nodes)) then
      call fix_enlist(enlist, nodes_translate)
    end if

    call allocate(mesh, dim, nnodes_local, nelements, loc)
    mesh%ndglno = reshape(enlist, (/loc * nelements/))
    deallocate(enlist)
    
    call allocate(positions, dim, mesh)
    positions%val = coords
    deallocate(coords)
    call deallocate(mesh)
    
    if(present(nodes)) then
      deallocate(nodes_translate)
    end if

  end function read_triangle_files

  subroutine fix_enlist(enlist, nodes_translate)
    integer, dimension(:, :), intent(inout) :: enlist
    integer, dimension(:), intent(in)       :: nodes_translate

    integer, dimension(size(enlist,1))      :: cell_nodes
    integer :: i, j, k

    do i = 1, size(enlist, 2)
      cell_nodes = enlist(:,i)
      do j = 1, size(nodes_translate)
        if (any(cell_nodes .eq. nodes_translate(j))) then
          do k = 1, size(enlist,1)
            if (cell_nodes(k) .eq. nodes_translate(j)) cell_nodes(k) = j
          end do
        end if
      end do
      enlist(:,i) = cell_nodes
    end do
  end subroutine fix_enlist
  
  subroutine read_halo_files(filename, nnodes, nodes)
    character(len = *), intent(in) :: filename
    integer, intent(out) :: nnodes
    integer, dimension(:), allocatable, intent(out) :: nodes

    call read_halo(trim(filename) // ".halo", nnodes, nodes)
  end subroutine read_halo_files

end module libsupermesh_read_triangle
