!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

module libsupermesh_halos_numbering

  use libsupermesh_data_structures
  use libsupermesh_fldebug
  use libsupermesh_futils
  use libsupermesh_halo_data_types
  use libsupermesh_halos_allocates
  use libsupermesh_halos_base
!  use halos_communications		! IAKOVOS commented out
!  use halos_debug			! IAKOVOS commented out
  use libsupermesh_mpi_interfaces
  use libsupermesh_parallel_tools
  use libsupermesh_quicksort

  implicit none
  
  private
  
!  public :: create_global_to_universal_numbering, &
!    & has_global_to_universal_numbering, universal_numbering_count, &
!    & halo_universal_number, halo_universal_numbers, get_universal_numbering, &
!    & get_universal_numbering_inverse, set_halo_universal_number, &
!    & ewrite_universal_numbers, valid_global_to_universal_numbering
  public :: create_global_to_universal_numbering, &
    & halo_universal_number, &
    & has_global_to_universal_numbering

  interface halo_universal_number
     module procedure halo_universal_number, halo_universal_number_vector
  end interface
    
!  interface get_universal_numbering
!     module procedure get_universal_numbering, get_universal_numbering_multiple_components
!  end interface

contains

  subroutine create_global_to_universal_numbering(halo, local_only)
    !!< Create the global to universal node numbering, and cache it on the halo
    !!< 
    !!< If local_only is present and .true. then only the universal numbers
    !!< for the owned nodes will be calculated. This is required when the
    !!< halos are not yet consistent and the universal numbers are to be
    !!< used to coordinate the halos.

    type(halo_type), intent(inout) :: halo
    logical, intent(in), optional :: local_only

    !assert(halo_valid_for_communication(halo))
    assert(.not. pending_communication(halo))

    if(has_global_to_universal_numbering(halo)) return

    select case(halo_ordering_scheme(halo))
      case(HALO_ORDER_GENERAL)
        call create_global_to_universal_numbering_order_general(halo, local_only)
      case(HALO_ORDER_TRAILING_RECEIVES)
        call create_global_to_universal_numbering_order_trailing_receives&
             (halo, local_only)
      case default
        FLAbort("Unrecognised halo ordering scheme")
    end select

#ifdef DDEBUG
    if(.not. present_and_true(local_only)) then
      assert(valid_global_to_universal_numbering(halo))
    end if
#endif

  end subroutine create_global_to_universal_numbering

  function has_global_to_universal_numbering(halo) result(has_gnn_to_unn)
    !!< Return whether the supplied halo has global to universal node numbering
    !!< data
    
    type(halo_type), intent(in) :: halo
    
    logical :: has_gnn_to_unn
   
    select case(halo_ordering_scheme(halo))
      case(HALO_ORDER_GENERAL)
        has_gnn_to_unn = associated(halo%owned_nodes_unn_base)
#ifdef DDEBUG
        if(has_gnn_to_unn) then
          assert(associated(halo%gnn_to_unn))
          assert(halo%my_owned_nodes_unn_base>=0)
        end if
#endif
      case(HALO_ORDER_TRAILING_RECEIVES)
        has_gnn_to_unn = associated(halo%receives_gnn_to_unn)
      case default
        FLAbort("Unrecognised halo ordering scheme")
    end select

  end function has_global_to_universal_numbering
  
  function halo_universal_number(halo, global_number) result(unn)
    !!< For the supplied halo, return the corresponding universal node number
    !!< for the supplied global node number
    
    type(halo_type), intent(in) :: halo
    integer, intent(in) :: global_number
    
    integer :: unn
    
    assert(global_number > 0)

    select case(halo_ordering_scheme(halo))
      case(HALO_ORDER_GENERAL)
        unn = halo_universal_number_order_general(halo, global_number) 
      case(HALO_ORDER_TRAILING_RECEIVES)
        unn = halo_universal_number_order_trailing_receives(halo, global_number) 
      case default
        FLAbort("Unrecognised halo ordering scheme")
    end select
    
  end function halo_universal_number

  function halo_universal_number_vector(halo, global_number) result(unn)
    !!< Version of halo_universal_number which returns a vector of
    !!< universal numbers corresponding to the supplied vector of global
    !!< numbers.    
    type(halo_type), intent(in) :: halo
    integer, dimension(:), intent(in) :: global_number

    integer, dimension(size(global_number)) :: unn
    
    integer :: i

    do i = 1, size(global_number)
       unn(i) = halo_universal_number(halo, global_number(i))
    end do

  end function halo_universal_number_vector
  
  function halo_universal_number_order_general(halo, global_number) result(unn)
    type(halo_type), intent(in) :: halo
    integer, intent(in) :: global_number
    
    integer :: unn

    assert(has_global_to_universal_numbering(halo))
    if (global_number<=size(halo%gnn_to_unn)) then
       unn = halo%gnn_to_unn(global_number)
    else
       unn = -1
    end if

  end function halo_universal_number_order_general
  
  function halo_universal_number_order_trailing_receives(halo, global_number) result(unn)
    type(halo_type), intent(in) :: halo
    integer, intent(in) :: global_number

    integer :: unn
    
    assert(has_global_to_universal_numbering(halo))

    
    if(global_number <= halo_nowned_nodes(halo)) then
      unn = halo%my_owned_nodes_unn_base + global_number
    else if(global_number - halo_nowned_nodes(halo) > size(halo%receives_gnn_to_unn)) then
      unn = - 1
    else
      unn = halo%receives_gnn_to_unn(global_number - halo_nowned_nodes(halo))
    end if
   
  end function halo_universal_number_order_trailing_receives

end module libsupermesh_halos_numbering
