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

module libsupermesh_halos_repair

  use libsupermesh_fields_data_types
  use libsupermesh_fields_base
  use libsupermesh_parallel_tools
  use libsupermesh_fldebug
  use libsupermesh_halo_data_types
  use libsupermesh_halos_base
!  use halos_debug		! IAKOVOS commented out
  use libsupermesh_halos_numbering
!  use halos_ownership		! IAKOVOS commented out
  use libsupermesh_mpi_interfaces
  use libsupermesh_quicksort

  implicit none

  private

!  public :: reorder_halo, reorder_l1_from_l2_halo, reorder_element_halo, &
!    & reorder_halo_receives, reorder_halo_from_element_halo
  public :: reorder_halo_from_element_halo
  
!  interface reorder_halo 
!    module procedure reorder_halo_vector, reorder_halo_halo
!  end interface reorder_halo 
  
contains

  subroutine reorder_halo_from_element_halo(node_halo, element_halo, mesh)
    !!< Using the order information in the element halo, reorder the sends
    !!< and receives in halo into a consistent order.
    !!<
    !!< This has the side effect of also defining the universal numbering on
    !!< node_halo.
    
    type(halo_type), intent(inout) :: node_halo
    type(halo_type), intent(in) :: element_halo
    type(mesh_type), intent(in) :: mesh

    integer :: p, nprocs, n
    integer, dimension(:), allocatable ::  global_numbers, order

! IAKOVOS commented out
    FLAbort("reorder_halo_from_element_halo: Code Commented out")
    
  end subroutine reorder_halo_from_element_halo
    
end module libsupermesh_halos_repair
