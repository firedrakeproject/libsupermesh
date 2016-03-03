#include "libsupermesh_debug.h"

module libsupermesh

  use libsupermesh_intersection_finder
  use libsupermesh_parallel_supermesh, only : parallel_supermesh
  use libsupermesh_supermesh

  implicit none

  public
  
  integer, parameter :: version_major = libsupermesh_VERSION_MAJOR
  integer, parameter :: version_minor = libsupermesh_VERSION_MINOR
  integer, parameter :: version_subminor = libsupermesh_VERSION_SUBMINOR
#ifdef libsupermesh_VERSION_RELEASE
  logical, parameter :: version_release = .true.
#else
  logical, parameter :: version_release = .false.
#endif
  character(len = *), parameter :: version = libsupermesh_VERSION_STRING

end module libsupermesh
