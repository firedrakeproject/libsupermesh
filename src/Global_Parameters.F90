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

module libsupermesh_global_parameters
  !!< This routine exists to save us all from argument list hell!
  !!<
  !!< All the global parameters which don't change while fluidity is running
  !!< should live here. I am building this up as I encounter more parameters
  !!< in the code. It would be great if others did the same.
  !!<
  !!< The correct syntax for accessing this module is:
  !!<
  !!< use global_parameters, only: parameter1, parameter2 ...
  !!<
  !!< Try to only use the parameters which are needed locally.

  implicit none
        
  !------------------------------------------------------------------------
  ! Precision parameters
  !------------------------------------------------------------------------
  !! Number of digits past the decimal point for a real
  integer, parameter :: real_digits_10 = precision(0.0)

  ! real variable declarations of the form:
  !   real*4 :: real_var
  ! are not portable. Use these instead:
  !   real(real_4) :: real_val
  integer, parameter ::real_8 = kind(0.0D0)

end module libsupermesh_global_parameters
