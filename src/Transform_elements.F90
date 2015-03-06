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
module libsupermesh_transform_elements
  ! Module to calculate element transformations from local to physical
  ! coordinates.
!  use quadrature			! IAKOVOS commented out
  use libsupermesh_elements
!  use vector_tools			! IAKOVOS commented out
  use libsupermesh_parallel_tools, only: abort_if_in_parallel_region
!  use fields_base			! IAKOVOS commented out
!  use cv_faces, only: cv_faces_type	! IAKOVOS commented out
!  use eventcounter			! IAKOVOS commented out
!  use memory_diagnostics		! IAKOVOS commented out
  
  implicit none
  
end module libsupermesh_transform_elements
