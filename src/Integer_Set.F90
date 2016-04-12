!  Copyright (C) 2016 The University of Edinburgh
!
!  The file is part of libsupermesh
!    
!  This library is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation;
!  version 2.1 of the License.
!
!  This library is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public
!  License along with this library; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

! The following code is derived from
! femtools/Integer_set.F90 in Fluidity git revision
! 4e6c1d2b022df3a519cdec120fad28e60d1b08d9 (dated 2015-02-25)
!
! Fluidity copyright information (note that AUTHORS mentioned in the following
! has been renamed to fluidity_AUTHORS):
!
!  Copyright (C) 2006 Imperial College London and others.
!  
!  Please see the AUTHORS file in the main source directory for a full list
!  of copyright holders.
!
!  Prof. C Pain
!  Applied Modelling and Computation Group
!  Department of Earth Science and Engineering
!  Imperial College London
!
!  amcgsoftware@imperial.ac.uk
!  
!  This library is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation,
!  version 2.1 of the License.
!
!  This library is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public
!  License along with this library; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!  USA

#include "libsupermesh_debug.h"

module libsupermesh_integer_set

  use iso_c_binding, only : c_ptr
  
  use libsupermesh_debug, only : abort_pinpoint

  implicit none
  
  private
  
  type integer_set
    type(c_ptr) :: address
  end type integer_set

  type integer_set_vector
     type(integer_set), dimension(:), pointer :: sets
  end type integer_set_vector

  interface
    subroutine integer_set_create_c(i) bind(c, name = "libsupermesh_integer_set_create")
      use iso_c_binding, only : c_ptr
      implicit none
      type(c_ptr) :: i
    end subroutine integer_set_create_c

    subroutine integer_set_delete_c(i) bind(c, name = "libsupermesh_integer_set_delete")
      use iso_c_binding, only : c_ptr
      implicit none
      type(c_ptr) :: i
    end subroutine integer_set_delete_c

    subroutine integer_set_insert_c(i, v, changed) bind(c, name = "libsupermesh_integer_set_insert")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int) :: v
      integer(kind = c_int) :: changed
    end subroutine integer_set_insert_c

    subroutine integer_set_length_c(i, l) bind(c, name = "libsupermesh_integer_set_length")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int) :: l
    end subroutine integer_set_length_c

    subroutine integer_set_fetch_c(i, idx, v) bind(c, name = "libsupermesh_integer_set_fetch")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int) :: idx
      integer(kind = c_int) :: v
    end subroutine integer_set_fetch_c

    subroutine integer_set_remove_c(i, v, status) bind(c, name = "libsupermesh_integer_set_remove")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int) :: v
      integer(kind = c_int) :: status
    end subroutine integer_set_remove_c

    subroutine integer_set_has_value_c(i, v, present) bind(c, name = "libsupermesh_integer_set_has_value")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int) :: v
      integer(kind = c_int) :: present
    end subroutine integer_set_has_value_c
  end interface

  interface allocate
    module procedure integer_set_allocate_single, integer_set_allocate_vector
  end interface allocate

  interface insert
    module procedure integer_set_insert, integer_set_insert_multiple, &
      integer_set_insert_set
  end interface insert

  interface deallocate
    module procedure integer_set_delete_single, integer_set_delete_vector
  end interface deallocate

  interface has_value
    module procedure integer_set_has_value, integer_set_has_value_multiple
  end interface has_value

  interface key_count
    module procedure integer_set_length_single, integer_set_length_vector
  end interface key_count

  interface fetch
    module procedure integer_set_fetch
  end interface fetch

  interface remove
    module procedure integer_set_remove
  end interface remove
  
  interface copy
    module procedure integer_set_copy, integer_set_copy_multiple
  end interface copy
  
  interface set_intersection
    module procedure set_intersection_two, set_intersection_multiple
  end interface set_intersection
  
  public :: integer_set, allocate, deallocate, has_value, key_count, fetch, &
    & insert, set_complement, set_intersection, set_minus, remove, copy, &
    & integer_set_vector

contains 
  
  subroutine integer_set_allocate_single(iset)
    type(integer_set), intent(out) :: iset
    iset = integer_set_create()
  end subroutine integer_set_allocate_single
  
  subroutine integer_set_allocate_vector(iset)
    type(integer_set), dimension(:), intent(out) :: iset
    
    integer :: i
    
    do i = 1, size(iset)
      call allocate(iset(i))
    end do
  
  end subroutine integer_set_allocate_vector

  function integer_set_create() result(iset)
    type(integer_set) :: iset
    call integer_set_create_c(iset%address)
  end function integer_set_create

  subroutine integer_set_delete_single(iset)
    type(integer_set), intent(inout) :: iset
    call integer_set_delete_c(iset%address)
  end subroutine integer_set_delete_single
  
  subroutine integer_set_delete_vector(iset)
    type(integer_set), dimension(:), intent(inout) :: iset
    
    integer :: i
    
    do i = 1, size(iset)
      call deallocate(iset(i))
    end do
    
  end subroutine integer_set_delete_vector

  subroutine integer_set_insert(iset, val, changed)
    type(integer_set), intent(inout) :: iset
    integer, intent(in) :: val
    logical, intent(out), optional :: changed
    integer :: lchanged

    call integer_set_insert_c(iset%address, val, lchanged)

    if (present(changed)) then
      changed = (lchanged == 1)
    end if
  end subroutine integer_set_insert

  subroutine integer_set_insert_multiple(iset, values)
    type(integer_set), intent(inout) :: iset
    integer, dimension(:), intent(in) :: values
    integer :: i

    do i=1,size(values)
      call insert(iset, values(i))
    end do
  end subroutine integer_set_insert_multiple

  subroutine integer_set_insert_set(iset, value_set)
    type(integer_set), intent(inout) :: iset
    type(integer_set), intent(in) :: value_set
    integer :: i

    do i=1, key_count(value_set)
      call insert(iset, fetch(value_set,i))
    end do
  end subroutine integer_set_insert_set
  
  function integer_set_length_single(iset) result(len)
    type(integer_set), intent(in) :: iset
    integer :: len

    call integer_set_length_c(iset%address, len)
  end function integer_set_length_single
  
  function integer_set_length_vector(iset) result(len)
    type(integer_set), dimension(:), intent(in) :: iset

    integer, dimension(size(iset)) :: len
    
    integer :: i
    
    do i = 1, size(iset)
      len(i) = key_count(iset(i))
    end do
  
  end function integer_set_length_vector

  function integer_set_fetch(iset, idx) result(val)
    type(integer_set), intent(in) :: iset
    integer, intent(in) :: idx
    integer :: val

    call integer_set_fetch_c(iset%address, idx, val)
  end function integer_set_fetch

  subroutine integer_set_remove(iset, val)
    type(integer_set), intent(in) :: iset
    integer, intent(in) :: val
    integer :: stat

    call integer_set_remove_c(iset%address, val, stat)
    assert(stat == 1)
  end subroutine integer_set_remove

  function integer_set_has_value(iset, val) result(bool)
    type(integer_set), intent(in) :: iset
    integer, intent(in) :: val
    logical :: bool

    integer :: lbool
    call integer_set_has_value_c(iset%address, val, lbool)
    bool = (lbool == 1)
  end function integer_set_has_value

  function integer_set_has_value_multiple(iset, val) result(bool)
    type(integer_set), intent(in) :: iset
    integer, dimension(:), intent(in) :: val
    logical, dimension(size(val)) :: bool
    
    integer:: i
    
    do i=1, size(val)
      bool(i)=integer_set_has_value(iset, val(i))
    end do
  end function integer_set_has_value_multiple
  
  subroutine set_complement(complement, universe, current)
    ! complement = universe \ current
    type(integer_set), intent(out) :: complement
    type(integer_set), intent(in) :: universe, current
    integer :: i, val

    call allocate(complement)
    do i=1,key_count(universe)
      val = fetch(universe, i)
      if (.not. has_value(current, val)) then
        call insert(complement, val)
      end if
    end do
  end subroutine set_complement

  subroutine set_intersection_two(intersection, A, B)
    ! intersection = A n B
    type(integer_set), intent(out) :: intersection
    type(integer_set), intent(in) :: A, B
    integer :: i, val

    call allocate(intersection)
    do i=1,key_count(A)
      val = fetch(A, i)
      if (has_value(B, val)) then
        call insert(intersection, val)
      end if
    end do
  end subroutine set_intersection_two

  subroutine set_intersection_multiple(intersection, isets)
    ! intersection = isets(i) n isets(j), forall i /= j
    type(integer_set), intent(out) :: intersection
    type(integer_set), dimension(:), intent(in) :: isets
    integer :: i
    
    type(integer_set) :: tmp_intersection, tmp_iset

    tmp_iset = isets(1)
    do i = 2, size(isets)
      call set_intersection(tmp_intersection, tmp_iset, isets(i))
      call copy(tmp_iset, tmp_intersection)
      call deallocate(tmp_intersection)
    end do
    intersection = tmp_iset
    
  end subroutine set_intersection_multiple
  
  subroutine integer_set_copy(iset_copy, iset)
    type(integer_set), intent(out) :: iset_copy
    type(integer_set), intent(in) :: iset
    
    integer :: i, val
    
    call allocate(iset_copy)
    
    do i = 1, key_count(iset)
      val = fetch(iset, i)
      call insert(iset_copy, val)
    end do
  
  end subroutine integer_set_copy

  subroutine integer_set_copy_multiple(iset_copy, iset)
    type(integer_set), dimension(:), intent(out) :: iset_copy
    type(integer_set), dimension(:), intent(in) :: iset
    
    integer :: n
    
    do n=1, size(iset)
      call copy(iset_copy(n), iset(n))
    end do
  
  end subroutine integer_set_copy_multiple
  
  subroutine set_minus(minus, A, B)
  ! minus = A \ B
    type(integer_set), intent(out) :: minus
    type(integer_set), intent(in) :: A, B
    integer :: i, val

    call allocate(minus)
    do i=1,key_count(A)
      val = fetch(A, i)
      if (.not. has_value(B, val)) then
        call insert(minus, val)
      end if
    end do
  end subroutine set_minus

end module libsupermesh_integer_set
