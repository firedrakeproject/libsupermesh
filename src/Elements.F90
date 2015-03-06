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
module libsupermesh_elements
  !!< This module provides derived types for finite elements and associated functions.
  use libsupermesh_element_numbering
!  use quadrature		! IAKOVOS commented out
  use libsupermesh_FLDebug
!  use polynomials		! IAKOVOS commented out
  use libsupermesh_reference_counting
  implicit none

  type element_type
     !!< Type to encode shape and quadrature information for an element.
     integer :: dim !! 2d or 3d?
     integer :: loc !! Number of nodes.
     integer :: ngi !! Number of gauss points.
     integer :: degree !! Polynomial degree of element.
     !! Shape functions: n is for the primitive function, dn is for partial derivatives, dn_s is for partial derivatives on surfaces. 
     !! n is loc x ngi, dn is loc x ngi x dim
     !! dn_s is loc x ngi x face x dim 
     real, pointer :: n(:,:)=>null(), dn(:,:,:)=>null()
     real, pointer :: n_s(:,:,:)=>null(), dn_s(:,:,:,:)=>null()
     !! Polynomials defining shape functions and their derivatives.
!     type(polynomial), dimension(:,:), pointer :: spoly=>null(), dspoly=>null()	! IAKOVOS commented out
     !! Link back to the node numbering used for this element.
     type(ele_numbering_type), pointer :: numbering=>null()
     !! Link back to the quadrature used for this element.
!     type(quadrature_type) :: quadrature				! IAKOVOS commented out
!     type(quadrature_type), pointer :: surface_quadrature=>null()	! IAKOVOS commented out
     !! Pointer to the superconvergence data for this element.
!     type(superconvergence_type), pointer :: superconvergence=>null()	! IAKOVOS commented out
     !! Pointer to constraints data for this element
!     type(constraints_type), pointer :: constraints=>null()		! IAKOVOS commented out
     !! Reference count to prevent memory leaks.
     type(refcount_type), pointer :: refcount=>null()
     !! Dummy name to satisfy reference counting
     character(len=0) :: name
  end type element_type
  
#include "Reference_count_interface_element_type.F90"
  
contains

#include "Reference_count_element_type.F90"

end module libsupermesh_elements
