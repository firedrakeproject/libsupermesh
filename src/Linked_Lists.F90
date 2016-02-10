#include "libsupermesh_debug.h"

module libsupermesh_linked_lists
  ! A module to provide linked lists and operations on them.
  use libsupermesh_debug
  implicit none

  ! Define a linked list for integers
  TYPE inode
     INTEGER :: value               
     TYPE (inode), POINTER :: next=>null() ! next node
  END TYPE inode
  
  TYPE ilist
     integer :: length=0
     TYPE (inode), POINTER :: firstnode=>null()
     type(inode), pointer ::  lastnode => null()
  END TYPE ilist
  
  interface insert_ascending
     module procedure iinsert_ascending
  end interface

  interface has_value
     module procedure ihas_value
  end interface

  interface deallocate
    module procedure flush_ilist, flush_ilist_v
  end interface

  interface insert
     module procedure iinsert
  end interface

  interface flush_list
     module procedure flush_ilist
  end interface
  
  interface flush_lists
    module procedure flush_ilist_array
  end interface flush_lists
  
  interface pop
     module procedure ipop
  end interface

  interface fetch
     module procedure ifetch
  end interface

  interface list2vector
     module procedure ilist2vector
  end interface

  interface pop_last
    module procedure ipop_last
  end interface

  interface size_intersection
    module procedure isize_intersection
  end interface

  interface has_value_sorted
    module procedure ihas_value_sorted
  end interface
  
  interface print_list
    module procedure iprint
  end interface
  
  interface intersect_ascending
    module procedure intersect_ascending_ilist
  end interface intersect_ascending

  interface copy
     module procedure copy_ilist, copy_ilist_array
  end interface

  interface maxval
    module procedure list_maxval
  end interface maxval

contains

  integer function list_maxval(list)
    type(ilist), intent(in) :: list
    type(inode), pointer :: node

    node => list%firstnode
    list_maxval = node%value
    do while (associated(node))
      list_maxval = max(list_maxval, node%value)
      node => node%next
    end do
  end function list_maxval

  logical function ihas_value(list, value)
    ! Check if the list contains the value.
    type(ilist), intent(in) :: list
    integer, intent(in) :: value
    
    type(inode), pointer :: node

    ihas_value = .false.

    node => list%firstnode
    do while (associated(node))
       if(value==node%value) then
          ihas_value = .true.
          return
       end if
       node => node%next
    end do
  end function ihas_value

  subroutine iinsert_ascending(list, value, discard)
    ! Insert value in list in such a position as to ensure that list remains
    ! in ascending order. This assumes that list is in ascending order.
    ! Duplicate values are discarded!
    type(ilist), intent(inout) :: list
    integer, intent(in) :: value
    
    logical, optional :: discard

    type(inode), pointer :: this_node, next_node
    integer :: pos

    ! Special case for zero length lists.
    if (list%length==0) then
       allocate(list%firstnode)

       list%firstnode%value=value
       ! The following should not be necessary
       list%firstnode%next=>null()
       
       list%length=1
       return
    end if
       
    this_node=>list%firstnode
    next_node=>list%firstnode%next

    ! Special case for a value smaller than the first value.
    if (value<list%firstnode%value) then
       allocate(list%firstnode)

       list%firstnode%next=>this_node
       
       list%firstnode%value=value
       
       list%length=list%length+1
       return
    end if

    ! initialise discard logical
    if (present(discard)) discard =.false.

    do pos=0,list%length
       if(this_node%value==value) then
          ! Discard duplicates.
          if (present(discard)) discard = .true.
          return
       end if

       if (.not.associated(next_node)) then
          ! We have hit then end of the chain. 
          allocate(this_node%next)

          if (this_node%value<value) then
             this_node%next%value=value
          else
             this_node%next%value=this_node%value
             this_node%value=value
          end if
          
          ! The following should not be necessary
          this_node%next%next=>null()
          
          list%length=list%length+1
          return
       end if

       ! Mid-chain. At this point we know this_node%value<value
       if (next_node%value>value) then
          ! Need to insert the value here.
          allocate(this_node%next)
          
          this_node%next%next=>next_node
       
          this_node%next%value=value
       
          list%length=list%length+1
          return
       end if

       ! Move along the chain.
       next_node=>next_node%next
       this_node=>this_node%next

    end do

    libsupermesh_abort("Walked off the end of the list. This can't happen.")

  end subroutine iinsert_ascending

  subroutine iinsert(list, i)
    type(ilist), intent(inout) :: list
    integer, intent(in) :: i
    type(inode), pointer :: node
    
    ! Special case for zero length lists.
    if (list%length==0) then
       allocate(list%firstnode)

       list%firstnode%value=i
       ! The following should not be necessary
       list%firstnode%next=>null()
       
       list%length=1
       list%lastnode => list%firstnode
       return
    end if
       
    node => list%lastnode
    allocate(node%next)
    node%next%value = i

    ! The following should not be necessary
    node%next%next => null()
          
    list%length = list%length+1
    list%lastnode => node%next
    return
  end subroutine iinsert

  subroutine flush_ilist(list)
    ! Remove all entries from a list.
    type(ilist), intent(inout) ::list

    integer :: i, tmp

    do i=1,list%length
       tmp=pop(list)
    end do

  end subroutine flush_ilist

  subroutine flush_ilist_v(lists)
    type(ilist), intent(inout), dimension(:) :: lists
    integer :: i

    do i=1,size(lists)
      call flush_ilist(lists(i))
    end do
  end subroutine flush_ilist_v
  
  subroutine flush_ilist_array(lists)
    ! Remove all entries from an array of lists
    
    type(ilist), dimension(:), intent(inout) :: lists
    
    integer :: i
    
    do i = 1, size(lists)
      call flush_list(lists(i))
    end do
  
  end subroutine flush_ilist_array

  function ipop(list)
    ! Pop the first value off list.
    integer :: ipop
    type(ilist), intent(inout) :: list
    
    type(inode), pointer :: firstnode
    
    ipop=list%firstnode%value
    
    firstnode=>list%firstnode
    
    list%firstnode=>firstnode%next

    deallocate(firstnode)

    list%length=list%length-1

  end function ipop

  function ipop_last(list)
    ! Pop the last value off list.
    integer :: ipop_last
    type(ilist), intent(inout) :: list
    
    type(inode), pointer :: prev_node => null(), node
    integer :: i

    node => list%firstnode
    do i=1,list%length-1
      prev_node => node
      node => node%next
    end do

    ipop_last = node%value
    deallocate(node)
    prev_node%next => null()
    list%length = list%length - 1
  end function ipop_last

  function ifetch(list, j)
    integer :: ifetch
    type(ilist), intent(inout) :: list
    integer, intent(in) :: j
    
    type(inode), pointer :: node
    integer :: i

    node => list%firstnode
    do i=1,j-1
      node => node%next
    end do

    ifetch = node%value
  end function ifetch

  function ilist2vector(list) result (vector)
    ! Return a vector containing the contents of ilist
    type(ilist), intent(in) :: list
    integer, dimension(list%length) :: vector
    
    type(inode), pointer :: this_node
    integer :: i

    this_node=>list%firstnode

    do i=1,list%length
       vector(i)=this_node%value
       
       this_node=>this_node%next
    end do

  end function ilist2vector

  function isize_intersection(listA, listB) result(x)
    type(ilist), intent(in) :: listA, listB
    type(inode), pointer :: nodeA, nodeB
    integer :: x

    x = 0
    nodeA => listA%firstnode
    do while(associated(nodeA))
      nodeB => listB%firstnode
      do while(associated(nodeB))
        if (nodeA%value == nodeB%value) then
          x = x + 1
          exit
        else
          nodeB => nodeB%next
        end if
      end do
      nodeA => nodeA%next
    end do
  end function isize_intersection

  function ihas_value_sorted(list, i) result(isin)
  ! This function assumes list is sorted
  ! in ascending order
    type(ilist), intent(in) :: list
    integer, intent(in) :: i
    type(inode), pointer :: node
    logical :: isin

    node => list%firstnode
    isin = .false.

    do while(associated(node))
      if (node%value > i) then
        return
      else if (node%value == i) then
        isin = .true.
        return
      end if
      node => node%next
    end do
  end function ihas_value_sorted

  subroutine iprint(list, priority)
    type(ilist), intent(in) :: list
    integer, intent(in) :: priority
    type(inode), pointer :: node

    ewrite(priority, *) "length: ", list%length

    node => list%firstnode
    do while (associated(node))
      ewrite(priority, *) " -- ", node%value
      node => node%next
    end do
  end subroutine
  
  function intersect_ascending_ilist(list1, list2) result(intersection)
    !!< Assumes that list1 and list2 are already sorted
    type(ilist), intent(in) :: list1
    type(ilist), intent(in) :: list2
    
    type(ilist) :: intersection
    
    type(inode), pointer :: node1 => null(), node2 => null()
    
    node1 => list1%firstnode
    node2 => list2%firstnode
    do while(associated(node1) .and. associated(node2))
      if(node1%value == node2%value) then
        call insert_ascending(intersection, node1%value)
        node1 => node1%next
        node2 => node2%next
      else
        if(node1%value < node2%value) then
          node1 => node1%next
        else
          node2 => node2%next
        end if
      end if
    end do
    
  end function intersect_ascending_ilist

  subroutine copy_ilist(copy_list, list)
    !!< Make a deep copy of list
    type(ilist), intent(out) :: copy_list
    type(ilist), intent(in) :: list
    
    type(inode), pointer :: node, copy_node

    if (list%length==0) return

    ! Special case the first entry
    node=>list%firstnode
    allocate(copy_list%firstnode)
    copy_list%firstnode%value=node%value
    copy_node=>copy_list%firstnode
    copy_list%length=1
    node=>node%next

    do while(associated(node))
       allocate(copy_node%next)
       copy_node=>copy_node%next
       copy_node%value=node%value
       
       copy_list%length=copy_list%length+1
       
       node=>node%next
    end do
    
  end subroutine copy_ilist

  subroutine copy_ilist_array(copy_lists, lists)
    !!< Make a deep copy of list
    type(ilist), dimension(:), intent(in) :: lists
    type(ilist), dimension(size(lists)), intent(out) :: copy_lists

    integer :: i

    do i=1,size(lists)
       call copy_ilist(copy_lists(i), lists(i))
    end do

  end subroutine copy_ilist_array

end module libsupermesh_linked_lists
