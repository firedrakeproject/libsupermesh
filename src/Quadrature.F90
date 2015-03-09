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
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
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
module libsupermesh_quadrature
  !!< This module implements quadrature of varying degrees for a number of
  !!< elements. Quadrature information is used to numerically evaluate
  !!< integrals over an element. 
  use libsupermesh_FLDebug
  use libsupermesh_reference_counting
!  use wandzura_quadrature		! IAKOVOS commented out
!  use grundmann_moeller_quadrature	! IAKOVOS commented out
!  use vector_tools			! IAKOVOS commented out
  implicit none

  private

  type quadrature_type
     !!< A data type which describes quadrature information. For most
     !!< developers, quadrature can be treated as an opaque data type which
     !!< will only be encountered when creating element_type variables to
     !!< represent shape functions.  
     integer :: dim !! Dimension of the elements for which quadrature
     !!< is required.  
     integer :: degree !! Degree of accuracy of quadrature. 
     integer :: vertices !! Number of vertices of the element.
     integer :: ngi !! Number of quadrature points.
     real, pointer :: weight(:)=>null() !! Quadrature weights.
     real, pointer :: l(:,:)=>null() !! Locations of quadrature points.
     character(len=0) :: name !! Fake name for reference counting.
     !! Reference count to prevent memory leaks.
     type(refcount_type), pointer :: refcount=>null()
     integer :: family
  end type quadrature_type

#include "Reference_count_interface_quadrature_type.F90"

  interface deallocate
     module procedure deallocate_quad
  end interface
  
  public deallocate, quadrature_type, &
       & operator(==), incref, addref, decref

contains

  subroutine deallocate_quad(quad,stat)
    !!< Since quadrature types contain pointers it is necessary to
    !!< explicitly deallocate them.
    !! The quadrature type to be deallocated.
    type(quadrature_type), intent(inout) :: quad
    !! Stat returns zero for successful completion and nonzero otherwise.
    integer, intent(out), optional :: stat

    integer :: lstat

    call decref(quad)
    if (has_references(quad)) then
       ! There are still references to this quad so we don't deallocate.
       return
    end if
    
    deallocate(quad%weight,quad%l, stat=lstat)
    
    if (present(stat)) then
       stat=lstat
    else if (lstat/=0) then
!       FLAbort("Error deallocating quad")		! ToDo
    end if

  end subroutine deallocate_quad
  
#include "Reference_count_quadrature_type.F90"
  
end module libsupermesh_quadrature
