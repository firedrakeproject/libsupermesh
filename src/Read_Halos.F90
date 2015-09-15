#include "fdebug.h"

module libsupermesh_read_halos

  use iso_c_binding

  use libsupermesh_fldebug

  implicit none
  
#include "mpif.h"
  
  private
  
  public :: halo_reader_reset, halo_reader_set_input, halo_reader_query_output, halo_reader_get_output
  public :: int_array, halo_type, read_halo, deallocate
  
  interface
    subroutine halo_reader_reset() bind(c, name = "cLibSuperMesh_halo_reader_reset")
      implicit none
    end subroutine halo_reader_reset
    
    function halo_reader_set_input(filename, filename_len, process, nprocs) bind(c, name = "cLibSuperMesh_halo_reader_set_input") result(errorCount)
      use iso_c_binding, only : c_char, c_int
      implicit none
      integer(kind = c_int) :: filename_len, process, nprocs
      character(kind = c_char) :: filename(filename_len)
      integer(kind = c_int) :: errorCount
    end function halo_reader_set_input
    
    subroutine halo_reader_query_output(level, nprocs, nsends, nreceives) bind(c, name = "cLibSuperMesh_halo_reader_query_output")
      use iso_c_binding, only : c_int
      implicit none
      integer(kind = c_int) :: level, nprocs
      integer(kind = c_int), dimension(nprocs) :: nsends, nreceives
    end subroutine halo_reader_query_output
    
    subroutine halo_reader_get_output(level, nprocs, nsends, nreceives, npnodes, send, recv) bind(c, name = "cLibSuperMesh_halo_reader_get_output")
      use iso_c_binding, only : c_int
      implicit none
      integer(kind = c_int) :: level, nprocs, npnodes
      integer(kind = c_int), dimension(nprocs) :: nsends, nreceives
      integer(kind = c_int), dimension(sum(nsends)) :: send      
      integer(kind = c_int), dimension(sum(nreceives)) :: recv
    end subroutine halo_reader_get_output
  end interface
  
  type int_array
    integer, dimension(:), pointer :: val
  end type
  
  type halo_type
    integer :: level
    integer :: process
    integer :: nprocs
    integer :: npnodes
    type(int_array), dimension(:), pointer :: send
    type(int_array), dimension(:), pointer :: recv
  end type halo_type
  
  interface deallocate
    module procedure deallocate_halo
  end interface deallocate
  
contains

  subroutine read_halo(filename, halo, level)
    character(len = *), intent(in) :: filename
    type(halo_type), intent(out) :: halo
    integer, optional, intent(in) :: level
    
    integer(kind = c_int) :: errorCount
    integer(kind = c_int), dimension(:), allocatable :: nsends, nreceives, send, recv
    integer :: i, ierr, index
    
    call MPI_Comm_rank(MPI_COMM_WORLD, halo%process, ierr)
    if(ierr /= MPI_SUCCESS) then
      FLAbort("Unable to determine process number")
    end if
    call MPI_Comm_size(MPI_COMM_WORLD, halo%nprocs, ierr)
    if(ierr /= MPI_SUCCESS) then
      FLAbort("Unable to determine number of processes")
    end if
    errorCount = halo_reader_set_input(trim(filename), len_trim(filename), halo%process, halo%nprocs)
    if(errorCount /= 0) then
      FLExit("Unable to read halo file '" // trim(filename) // "'")
    end if
    
    if(present(level)) then
      halo%level = level
    else
      halo%level = 2
    end if
    allocate(nsends(halo%nprocs), nreceives(halo%nprocs))
    call halo_reader_query_output(halo%level, halo%nprocs, nsends, nreceives)
    allocate(halo%send(halo%nprocs), halo%recv(halo%nprocs))
    do i = 1, halo%nprocs
      allocate(halo%send(i)%val(nsends(i)), halo%recv(i)%val(nreceives(i)))
    end do
    
    allocate(send(sum(nsends)), recv(sum(nreceives)))
    call halo_reader_get_output(halo%level, halo%nprocs, nsends, nreceives, halo%npnodes, send, recv)
    index = 1
    do i = 1, halo%nprocs
      halo%send(i)%val = send(index:index + nsends(i) - 1)
      index = index + nsends(i)
    end do
    index = 1
    do i = 1, halo%nprocs
      halo%recv(i)%val = recv(index:index + nreceives(i) - 1)
      index = index + nreceives(i)
    end do
    deallocate(nsends, nreceives, send, recv)
    
  end subroutine read_halo
  
  subroutine deallocate_halo(halo)
    type(halo_type), intent(inout) :: halo
    
    integer :: i
    
    do i = 1, halo%nprocs
      deallocate(halo%send(i)%val, halo%recv(i)%val)
    end do
    deallocate(halo%send, halo%recv)
    
  end subroutine deallocate_halo
  
end module libsupermesh_read_halos
