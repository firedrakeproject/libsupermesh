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

module libsupermesh_parallel_tools

  use libsupermesh_fldebug
  use libsupermesh_mpi_interfaces
  use libsupermesh_global_parameters, only: is_active_process, no_active_processes
  use iso_c_binding
#ifdef _OPENMP
  use omp_lib
#endif
  implicit none

  private

  public :: abort_if_in_parallel_region
  public :: getnprocs, getprocno, &
       isparallel, &
       MPI_COMM_FEMTOOLS_LIB
  
  integer(c_int), bind(c) :: MPI_COMM_FEMTOOLS_LIB = MPI_COMM_WORLD
  
contains

  function getprocno(communicator) result(procno)
    !!< This is a convenience routine which returns the MPI rank
    !!< number + 1 when MPI is being used and 1 otherwise.
  
    integer, optional, intent(in) :: communicator

    integer :: procno
#ifdef HAVE_MPI
    integer :: ierr, lcommunicator
    logical :: initialized
    
    call MPI_Initialized(initialized, ierr)
    if(initialized) then
       if(present(communicator)) then
          lcommunicator = communicator
       else
          lcommunicator = MPI_COMM_FEMTOOLS_LIB
       end if
       
       assert(valid_communicator(lcommunicator))
       call MPI_Comm_Rank(lcommunicator, procno, ierr)
       assert(ierr == MPI_SUCCESS)
       procno = procno + 1
    else
       procno = 1
    end if
#else
    procno = 1
#endif

  end function getprocno

  ! Abort run if we're in an OMP parallel region
  ! Call this routine at the start of functions that are known not to
  ! be thread safe (for example, populating caches) and should
  ! therefore never be called in a parallel region due to race
  ! conditions.
  subroutine abort_if_in_parallel_region()
#ifdef _OPENMP
    if (omp_in_parallel()) then
       FLAbort("Calling non-thread-safe code in OMP parallel region")
    endif
#else
    return
#endif
  end subroutine abort_if_in_parallel_region
  
  function getnprocs(communicator) result(nprocs)
    !!< This is a convience routine which returns the number of processes
    !!< in a communicator (default MPI_COMM_FEMTOOLS_LIB) when MPI is being used and 1
    !!< otherwise.
    
    integer, optional, intent(in) :: communicator

    integer :: nprocs
    
#ifdef HAVE_MPI
    integer :: ierr, lcommunicator
    logical :: initialized

    call MPI_Initialized(initialized, ierr)
    if(initialized) then
       if(present(communicator)) then
          assert(valid_communicator(communicator))
          lcommunicator = communicator
       else
          lcommunicator = MPI_COMM_FEMTOOLS_LIB
       end if
       
       assert(valid_communicator(lcommunicator))
       call MPI_Comm_Size(lcommunicator, nprocs, ierr)
       assert(ierr == MPI_SUCCESS)
    else
       nprocs = 1
    end if
#else
    nprocs = 1
#endif

  end function getnprocs
  
  logical function isparallel()
    !!< Return true if we are running in parallel, and false otherwise.
  
    isparallel = (getnprocs()>1)
    
  end function isparallel

end module libsupermesh_parallel_tools
