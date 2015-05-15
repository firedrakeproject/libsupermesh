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
  public :: alland, allmax, allsum, allmean, &
       getnprocs, getpinteger, getpreal, getprocno, getrank, &
       isparallel, parallel_filename, parallel_filename_len, &
       pending_communication, valid_communicator, next_mpi_tag, &
       MPI_COMM_FEMTOOLS_LIB
  
  integer(c_int), bind(c) :: MPI_COMM_FEMTOOLS_LIB = MPI_COMM_WORLD
  
  interface allmax
    module procedure allmax_integer, allmax_real
  end interface allmax
  
  interface allsum
    module procedure allsum_integer, allsum_real, allsum_integer_vector, &
      & allsum_real_vector
  end interface allsum

  interface parallel_filename_len
    module procedure parallel_filename_no_extension_len, &
      &  parallel_filename_with_extension_len
  end interface
  
  interface parallel_filename
    module procedure parallel_filename_no_extension, &
      & parallel_filename_with_extension
  end interface
  
  interface pending_communication
    module procedure pending_communication_communicator
  end interface pending_communication
  
contains

  integer function next_mpi_tag()
#ifdef HAVE_MPI
    integer, save::last_tag=0, tag_ub=0
    integer flag, ierr
    if(tag_ub==0) then
       call MPI_Attr_get(MPI_COMM_FEMTOOLS_LIB, MPI_TAG_UB, tag_ub, flag, ierr)
    end if

    last_tag = mod(last_tag+1, tag_ub)
    if(last_tag==0) then
       last_tag = last_tag+1
    end if
    next_mpi_tag = last_tag
#else
    next_mpi_tag = 1
#endif
  end function next_mpi_tag
  
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
  
  function getrank(communicator) result(rank)
    !!< This is a convience routine which returns the MPI rank
    !!< number of the process when MPI is being used and 0 otherwise.
  
    integer, optional, intent(in) :: communicator

    integer::rank
#ifdef HAVE_MPI
    integer :: ierr, lcommunicator

    if(present(communicator)) then
      lcommunicator = communicator
    else
      lcommunicator = MPI_COMM_FEMTOOLS_LIB
    end if

    assert(valid_communicator(lcommunicator))
    call MPI_Comm_Rank(lcommunicator, rank, ierr)
    assert(ierr == MPI_SUCCESS)
#else
    rank = 0
#endif

  end function getrank

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
  
  function getpinteger() result(pinteger)
    !!< This is a convience routine which returns the MPI integer type 
    !!< being used. If MPI is not being used Pinteger is set to -1
  
    integer :: pinteger
#ifdef HAVE_MPI
    pinteger = MPI_INTEGER
#else
    pinteger = -1
#endif

  end function getpinteger
  
  function getpreal() result(preal)
    !!< This is a convience routine which returns the MPI real type 
    !!< being used. If MPI is not being used PREAL is set to -1
  
    integer :: preal
#ifdef HAVE_MPI

#ifdef DOUBLEP
    preal = MPI_DOUBLE_PRECISION
#else
    preal = MPI_REAL
#endif

#else
    preal = -1
#endif

  end function getpreal
  
  function pending_communication_communicator(communicator) result(pending)
    !!< Return whether there is a pending communication for the supplied communicator.
    
    integer, optional, intent(in) :: communicator

    logical :: pending

    integer :: lcommunicator

#ifdef HAVE_MPI
    integer :: ierr, ipending
    
    if(present(communicator)) then
      lcommunicator = communicator
    else
      lcommunicator = MPI_COMM_FEMTOOLS_LIB
    end if    

    call mpi_iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, lcommunicator, ipending, MPI_STATUS_IGNORE, ierr)
    assert(ierr == MPI_SUCCESS)
    
    pending = (ipending /= 0)

    ! Note - removing this mpi_barrier could result in a false
    ! positive on another process.
    call mpi_barrier(lcommunicator, ierr)
    assert(ierr == MPI_SUCCESS)
#else
    pending = .false.
#endif
    
  end function pending_communication_communicator
  
  function valid_communicator(communicator) result(valid)
    !!< Return whether the supplied MPI communicator is valid
    
    integer, intent(in) :: communicator
    
    logical :: valid
    
#ifdef HAVE_MPI
    integer :: ierr, size

    call mpi_comm_size(communicator, size, ierr)

    valid = (ierr == MPI_SUCCESS)
#else
    valid = .false.
#endif

  end function valid_communicator
  
  subroutine alland(value, communicator)
    !!< And the logical value across all processes
    
    logical, intent(inout) :: value
    integer, optional, intent(in) :: communicator
    
#ifdef HAVE_MPI
    integer :: ierr, lcommunicator
    logical :: and
    
    if(present(communicator)) then
      lcommunicator = communicator
    else
      lcommunicator = MPI_COMM_FEMTOOLS_LIB
    end if
    
    if(isparallel()) then    
      assert(valid_communicator(lcommunicator))
      call mpi_allreduce(value, and, 1, MPI_LOGICAL, MPI_LAND, lcommunicator, ierr)
      assert(ierr == MPI_SUCCESS)
      value = and
    end if
#endif

  end subroutine alland
  
  subroutine allmax_integer(value, communicator)
    !!< Find the maxmimum value across all processes
  
    integer, intent(inout) :: value
    integer, optional, intent(in) :: communicator
    
#ifdef HAVE_MPI
    integer :: ierr, lcommunicator, maximum
    
    if(present(communicator)) then
      lcommunicator = communicator
    else
      lcommunicator = MPI_COMM_FEMTOOLS_LIB
    end if
    
    if(isparallel()) then
       assert(valid_communicator(lcommunicator))
       call MPI_Allreduce(value, maximum, 1, getpinteger(), MPI_MAX, lcommunicator, ierr)
       assert(ierr == MPI_SUCCESS)
       value = maximum
    end if
#endif
  
  end subroutine allmax_integer

  subroutine allmax_real(value, communicator)
    !!< Find the maxmimum value across all processes
  
    real, intent(inout) :: value
    integer, optional, intent(in) :: communicator
    
#ifdef HAVE_MPI
    integer :: ierr, lcommunicator
    real :: maximum
    
    if(present(communicator)) then
      lcommunicator = communicator
    else
      lcommunicator = MPI_COMM_FEMTOOLS_LIB
    end if
    
    if(isparallel()) then
       assert(valid_communicator(lcommunicator))
       call mpi_allreduce(value, maximum, 1, getpreal(), MPI_MAX, lcommunicator, ierr)
       assert(ierr == MPI_SUCCESS)
       value = maximum
    end if
#endif

  end subroutine allmax_real
  
  subroutine allsum_integer(value, communicator)
    !!< Sum the integer value across all processes
  
    integer, intent(inout) :: value
    integer, optional, intent(in) :: communicator
    
#ifdef HAVE_MPI
    integer :: ierr, lcommunicator
    integer :: sum

    if(present(communicator)) then
      lcommunicator = communicator
    else
      lcommunicator = MPI_COMM_FEMTOOLS_LIB
    end if

    if(isparallel()) then
       assert(valid_communicator(lcommunicator))
       sum = 0.0
       call MPI_Allreduce(value, sum, 1, getpinteger(), MPI_SUM, lcommunicator, ierr)
       assert(ierr == MPI_SUCCESS)
       value = sum
    end if
#endif

  end subroutine allsum_integer

  subroutine allsum_real(value, communicator)
    !!< Sum the real value across all processes
  
    real, intent(inout) :: value
    integer, optional, intent(in) :: communicator
    
#ifdef HAVE_MPI
    integer :: ierr, lcommunicator
    real :: sum

    if(present(communicator)) then
      lcommunicator = communicator
    else
      lcommunicator = MPI_COMM_FEMTOOLS_LIB
    end if

    if(isparallel()) then
       assert(valid_communicator(lcommunicator))
       sum = 0.0
       call MPI_Allreduce(value, sum, 1, getpreal(), MPI_SUM, lcommunicator, ierr)
       assert(ierr == MPI_SUCCESS)
       value = sum
    end if
#endif

  end subroutine allsum_real
  
  subroutine allmean(value, communicator)
    !!< Sum the real value across all processes
  
    real, intent(inout) :: value
    integer, optional, intent(in) :: communicator
    
    call allsum(value, communicator = communicator)
    value = value / getnprocs(communicator = communicator)
    
  end subroutine allmean

  subroutine allsum_integer_vector(value, communicator)
    !!< Sum the value across all processes
  
    integer, intent(inout) :: value(:)
    integer, optional, intent(in) :: communicator
    
#ifdef HAVE_MPI
    integer :: lcommunicator, ierr
    integer, dimension(size(value)) :: sum
    
    if(present(communicator)) then
      lcommunicator = communicator
    else
      lcommunicator = MPI_COMM_FEMTOOLS_LIB
    end if
    
    if(isparallel()) then
       assert(valid_communicator(lcommunicator))
       sum = 0
       call MPI_Allreduce(value, sum, size(value), getpinteger(), MPI_SUM, lcommunicator, ierr)
       assert(ierr == MPI_SUCCESS)
       value = sum
    end if
#endif

  end subroutine allsum_integer_vector

  subroutine allsum_real_vector(value, communicator)
    !!< Sum the value across all processes
  
    real, intent(inout) :: value(:)
    integer, optional, intent(in) :: communicator
    
#ifdef HAVE_MPI
    integer :: lcommunicator, ierr
    real, dimension(size(value)) :: sum
    
    if(present(communicator)) then
      lcommunicator = communicator
    else
      lcommunicator = MPI_COMM_FEMTOOLS_LIB
    end if
    
    if(isparallel()) then
       assert(valid_communicator(lcommunicator))
       sum = 0.0
       call MPI_Allreduce(value, sum, size(value), getpreal(), MPI_SUM, lcommunicator, ierr)
       assert(ierr == MPI_SUCCESS)
       value = sum
    end if
#endif

  end subroutine allsum_real_vector
  
  pure function parallel_filename_no_extension_len(filename) result(length)
    !!< Return the (maximum) length of a string containing:
    !!<   [filename]_[process number]

    character(len = *), intent(in) :: filename
   
    integer :: length
    
    length = len_trim(filename) + 1 + floor(log10(real(huge(0)))) + 1
  
  end function parallel_filename_no_extension_len
  
  function parallel_filename_no_extension(filename) result(pfilename)
    !!< Return a string containing:
    !!<   [filename]-[process-number]
    !!< Note that is it important to trim the returned string.

    character(len = *), intent(in) :: filename
    
    character(len = parallel_filename_len(filename)) :: pfilename

    if (is_active_process .and. no_active_processes == 1) then
      write(pfilename, "(a)") trim(filename)
    else
      write(pfilename, "(a, i0)") trim(filename) // "_", getrank()
    end if
      
  end function parallel_filename_no_extension
  
  pure function parallel_filename_with_extension_len(filename, extension) result(length)
    !!< Return the (maximum) length of a string containing:
    !!<   [filename]_[process number].[extension]

    character(len = *), intent(in) :: filename
    character(len = *), intent(in) :: extension

    integer :: length

    length = parallel_filename_len(filename) + len_trim(extension)

  end function parallel_filename_with_extension_len
  
  function parallel_filename_with_extension(filename, extension)  result(pfilename)
    !!< Return a string containing:
    !!<   [filename]-[process-number][extension]
    !!< Note that is it important to trim the returned string.

    character(len = *), intent(in) :: filename
    character(len = *), intent(in) :: extension

    character(len = parallel_filename_len(filename, extension)) :: pfilename

    pfilename = trim(parallel_filename(filename)) // trim(extension)

  end function parallel_filename_with_extension

end module libsupermesh_parallel_tools
