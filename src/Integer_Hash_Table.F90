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
! femtools/Integer_hash_table.F90 in Fluidity git revision
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

module libsupermesh_integer_hash_table

  use iso_c_binding, only : c_int, c_ptr
  
  implicit none

  private
  
  type integer_hash_table
    type(c_ptr), pointer :: ptr
  end type integer_hash_table

  interface
    subroutine cinteger_hash_table_new(i) bind(c, name = "libsupermesh_integer_hash_table_new")
      use iso_c_binding, only : c_ptr
      implicit none
      type(c_ptr) :: i
    end subroutine cinteger_hash_table_new

    subroutine cinteger_hash_table_delete(i) bind(c, name = "libsupermesh_integer_hash_table_delete")
      use iso_c_binding, only : c_ptr
      implicit none
      type(c_ptr) :: i
    end subroutine cinteger_hash_table_delete

    subroutine cinteger_hash_table_insert(i, key, value) bind(c, name = "libsupermesh_integer_hash_table_insert")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int), value :: key
      integer(kind = c_int), value :: value
    end subroutine cinteger_hash_table_insert

    subroutine cinteger_hash_table_size(i, size) bind(c, name = "libsupermesh_integer_hash_table_size")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int) :: size
    end subroutine cinteger_hash_table_size

    subroutine cinteger_hash_table_fetch(i, key, value) bind(c, name = "libsupermesh_integer_hash_table_fetch")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int), value :: key
      integer(kind = c_int) :: value
    end subroutine cinteger_hash_table_fetch

    subroutine cinteger_hash_table_remove(i, key) bind(c, name = "libsupermesh_integer_hash_table_remove")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int), value :: key
    end subroutine cinteger_hash_table_remove

    subroutine cinteger_hash_table_has_key(i, key, present) bind(c, name = "libsupermesh_integer_hash_table_has_key")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int), value :: key
      integer(kind = c_int) :: present
    end subroutine cinteger_hash_table_has_key

    subroutine cinteger_hash_table_fetch_pair(i, index, key, value) bind(c, name = "libsupermesh_integer_hash_table_fetch_pair")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int), value :: index
      integer(kind = c_int) :: key
      integer(kind = c_int) :: value
    end subroutine cinteger_hash_table_fetch_pair
  end interface

  interface allocate
    module procedure allocate_integer_hash_table
  end interface allocate
  
  interface deallocate
    module procedure deallocate_integer_hash_table
  end interface deallocate

  interface insert
    module procedure integer_hash_table_insert
  end interface insert

  interface key_count
    module procedure integer_hash_table_size
  end interface key_count

  interface fetch
    module procedure integer_hash_table_fetch, integer_hash_table_fetch_rank_1
  end interface fetch

  interface remove
    module procedure integer_hash_table_remove
  end interface remove

  interface has_key
    module procedure integer_hash_table_has_key
  end interface has_key

  interface fetch_pair
    module procedure integer_hash_table_fetch_pair
  end interface fetch_pair

  public :: integer_hash_table, allocate, deallocate, insert, key_count, &
    & fetch, remove, has_key, fetch_pair

contains 

  subroutine allocate_integer_hash_table(ihash)
    type(integer_hash_table), intent(out) :: ihash

    allocate(ihash%ptr)
    call cinteger_hash_table_new(ihash%ptr)

  end subroutine allocate_integer_hash_table

  subroutine deallocate_integer_hash_table(ihash)
    type(integer_hash_table), intent(inout) :: ihash

    call cinteger_hash_table_delete(ihash%ptr)
    deallocate(ihash%ptr)

  end subroutine deallocate_integer_hash_table

  subroutine integer_hash_table_insert(ihash, key, value)
    type(integer_hash_table), intent(inout) :: ihash
    integer, intent(in) :: key
    integer, intent(in) :: value

    call cinteger_hash_table_insert(ihash%ptr, key, value)

  end subroutine integer_hash_table_insert

  function integer_hash_table_size(ihash) result(s)
    type(integer_hash_table), intent(inout) :: ihash

    integer :: s

    call cinteger_hash_table_size(ihash%ptr, s)

  end function integer_hash_table_size

  function integer_hash_table_fetch(ihash, key) result(value)
    type(integer_hash_table), intent(inout) :: ihash
    integer, intent(in) :: key

    integer :: value

    call cinteger_hash_table_fetch(ihash%ptr, key, value)

  end function integer_hash_table_fetch
  
  function integer_hash_table_fetch_rank_1(ihash, keys) result(values)
    type(integer_hash_table), intent(inout) :: ihash
    integer, dimension(:), intent(in) :: keys

    integer, dimension(size(keys)) :: values

    integer :: i

    do i = 1, size(keys)
      call cinteger_hash_table_fetch(ihash%ptr, keys(i), values(i))
    end do

  end function integer_hash_table_fetch_rank_1

  subroutine integer_hash_table_remove(ihash, key)
    type(integer_hash_table), intent(inout) :: ihash
    integer, intent(in) :: key

    call cinteger_hash_table_remove(ihash%ptr, key)

  end subroutine integer_hash_table_remove

  function integer_hash_table_has_key(ihash, key) result(present)
    type(integer_hash_table), intent(inout) :: ihash
    integer, intent(in) :: key

    logical :: present

    integer(kind = c_int) :: lpresent

    call cinteger_hash_table_has_key(ihash%ptr, key, lpresent)
    present = (lpresent /= 0)

  end function integer_hash_table_has_key

  subroutine integer_hash_table_fetch_pair(ihash, index, key, value)
    type(integer_hash_table), intent(inout) :: ihash
    integer, intent(in) :: index
    integer, intent(out) :: key
    integer, intent(out) :: value

    call cinteger_hash_table_fetch_pair(ihash%ptr, index, key, value)

  end subroutine integer_hash_table_fetch_pair

end module libsupermesh_integer_hash_table
