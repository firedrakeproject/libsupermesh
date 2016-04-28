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

  use iso_c_binding, only : c_ptr
  
  use libsupermesh_debug, only : abort_pinpoint

  implicit none

  private
  
  type integer_hash_table
    type(c_ptr) :: address
  end type integer_hash_table

  interface
    subroutine integer_hash_table_create_c(i) bind(c, name = "libsupermesh_integer_hash_table_create")
      use iso_c_binding, only : c_ptr
      implicit none
      type(c_ptr) :: i
    end subroutine integer_hash_table_create_c

    subroutine integer_hash_table_delete_c(i) bind(c, name = "libsupermesh_integer_hash_table_delete")
      use iso_c_binding, only : c_ptr
      implicit none
      type(c_ptr) :: i
    end subroutine integer_hash_table_delete_c

    subroutine integer_hash_table_insert_c(i, k, v) bind(c, name = "libsupermesh_integer_hash_table_insert")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int) :: k
      integer(kind = c_int) :: v
    end subroutine integer_hash_table_insert_c

    subroutine integer_hash_table_length_c(i, l) bind(c, name = "libsupermesh_integer_hash_table_length")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int) :: l
    end subroutine integer_hash_table_length_c

    subroutine integer_hash_table_fetch_c(i, k, v) bind(c, name = "libsupermesh_integer_hash_table_fetch")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int) :: k
      integer(kind = c_int) :: v
    end subroutine integer_hash_table_fetch_c

    subroutine integer_hash_table_remove_c(i, k, status) bind(c, name = "libsupermesh_integer_hash_table_remove")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int) :: k
      integer(kind = c_int) :: status
    end subroutine integer_hash_table_remove_c

    subroutine integer_hash_table_has_key_c(i, k, present) bind(c, name = "libsupermesh_integer_hash_table_has_key")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int) :: k
      integer(kind = c_int) :: present
    end subroutine integer_hash_table_has_key_c

    subroutine integer_hash_table_fetch_pair_c(i, idx, k, v) bind(c, name = "libsupermesh_integer_hash_table_fetch_pair")
      use iso_c_binding, only : c_int, c_ptr
      implicit none
      type(c_ptr) :: i
      integer(kind = c_int) :: idx
      integer(kind = c_int) :: k
      integer(kind = c_int) :: v
    end subroutine integer_hash_table_fetch_pair_c
  end interface

  interface allocate
    module procedure integer_hash_table_allocate
  end interface allocate

  interface insert
    module procedure integer_hash_table_insert
  end interface insert

  interface remove
    module procedure integer_hash_table_remove
  end interface remove

  interface deallocate
    module procedure integer_hash_table_delete
  end interface deallocate

  interface has_key
    module procedure integer_hash_table_has_key
  end interface has_key

  interface key_count
    module procedure integer_hash_table_length
  end interface key_count

  interface fetch
    module procedure integer_hash_table_fetch, integer_hash_table_fetch_v
  end interface fetch

  interface fetch_pair
    module procedure integer_hash_table_fetch_pair
  end interface fetch_pair

  interface copy
    module procedure integer_hash_table_copy
  end interface copy

  public :: integer_hash_table, allocate, deallocate, has_key, key_count, &
    & fetch, insert, fetch_pair, remove, copy

contains 

  subroutine integer_hash_table_copy(ihash_copy, ihash)
    type(integer_hash_table), intent(out) :: ihash_copy
    type(integer_hash_table), intent(in) :: ihash

    integer :: ind, key, key_val

    call allocate(ihash_copy)
    do ind = 1, key_count(ihash)
      call fetch_pair(ihash, ind, key, key_val)
      call insert(ihash_copy, key, key_val)
    end do

  end subroutine integer_hash_table_copy

  subroutine integer_hash_table_allocate(ihash)
    type(integer_hash_table), intent(out) :: ihash
    ihash = integer_hash_table_create()
  end subroutine integer_hash_table_allocate

  function integer_hash_table_create() result(ihash)
    type(integer_hash_table) :: ihash
    call integer_hash_table_create_c(ihash%address)
  end function integer_hash_table_create

  subroutine integer_hash_table_delete(ihash)
    type(integer_hash_table), intent(inout) :: ihash
    call integer_hash_table_delete_c(ihash%address)
  end subroutine integer_hash_table_delete

  subroutine integer_hash_table_insert(ihash, key, val)
    type(integer_hash_table), intent(inout) :: ihash
    integer, intent(in) :: key, val

    call integer_hash_table_insert_c(ihash%address, key, val)
  end subroutine integer_hash_table_insert

  function integer_hash_table_length(ihash) result(len)
    type(integer_hash_table), intent(in) :: ihash
    integer :: len

    call integer_hash_table_length_c(ihash%address, len)
  end function integer_hash_table_length

  function integer_hash_table_fetch(ihash, key) result(val)
    type(integer_hash_table), intent(in) :: ihash
    integer, intent(in) :: key
    integer :: val

    call integer_hash_table_fetch_c(ihash%address, key, val)
  end function integer_hash_table_fetch

  subroutine integer_hash_table_remove(ihash, key)
    type(integer_hash_table), intent(inout) :: ihash
    integer, intent(in) :: key
    integer :: stat

    call integer_hash_table_remove_c(ihash%address, key, stat)
    assert(stat == 1)
  end subroutine integer_hash_table_remove

  function integer_hash_table_fetch_v(ihash, keys) result(vals)
    type(integer_hash_table), intent(in) :: ihash
    integer, intent(in), dimension(:) :: keys
    integer, dimension(size(keys)) :: vals
    integer :: i

    do i=1,size(keys)
      call integer_hash_table_fetch_c(ihash%address, keys(i), vals(i))
    end do
  end function integer_hash_table_fetch_v

  function integer_hash_table_has_key(ihash, key) result(bool)
    type(integer_hash_table), intent(in) :: ihash
    integer, intent(in) :: key
    logical :: bool

    integer :: lbool
    call integer_hash_table_has_key_c(ihash%address, key, lbool)
    bool = (lbool == 1)
  end function integer_hash_table_has_key

  subroutine integer_hash_table_fetch_pair(ihash, idx, key, val)
    type(integer_hash_table), intent(in) :: ihash
    integer, intent(in) :: idx
    integer, intent(out) :: key, val

    call integer_hash_table_fetch_pair_c(ihash%address, idx, key, val)
  end subroutine integer_hash_table_fetch_pair
  
end module libsupermesh_integer_hash_table
