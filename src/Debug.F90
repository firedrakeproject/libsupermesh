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
! debug/Debug.F90 in Fluidity git revision
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

module libsupermesh_debug
      
  use libsupermesh_debug_parameters, only : current_debug_level, &
    & debug_error_unit, debug_log_unit
  
  implicit none
  
#include <mpif.h>
  
  private
  
  public :: current_debug_level, debug_unit, abort_pinpoint

  interface print_backtrace
    subroutine libsupermesh_print_backtrace(max_size) bind(c)
      use iso_c_binding, only : c_int
      implicit none
      integer(kind = c_int) :: max_size
    end subroutine libsupermesh_print_backtrace
  end interface print_backtrace
  
contains
  
  function debug_unit(priority)    
    integer, intent(in) :: priority
    
    integer :: debug_unit
    
    if(priority < 1) then
       debug_unit = debug_error_unit
    else
       debug_unit = debug_log_unit
    end if
    
  end function debug_unit

  subroutine abort_pinpoint(error, file, line_number)
    character(len = *), intent(in) :: error
    character(len = *), intent(in) :: file
    integer, intent(in) :: line_number
    
    integer :: ierr, mpi_init
    
    call MPI_Initialized(mpi_init, ierr)
    ewrite(-1, "(a)")      "*** libsupermesh error ***"
    ewrite(-1, "(a,i0,a)") "Source location: (" // trim(file) // ",", line_number, ")"
    ewrite(-1, "(a)")      "Error message: " // trim(error)
    call print_backtrace(max_size = 64)
    if(mpi_init /= 0) call MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER, ierr)
    stop 1
    
  end subroutine abort_pinpoint

end module libsupermesh_debug
