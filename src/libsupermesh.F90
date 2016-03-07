#include "libsupermesh_debug.h"

module libsupermesh

  use libsupermesh_intersection_finder
  use libsupermesh_parallel_supermesh, only : parallel_supermesh
  use libsupermesh_supermesh

  implicit none

  public
  
  integer, parameter :: version_major = LIBSUPERMESH_VERSION_MAJOR
  integer, parameter :: version_minor = LIBSUPERMESH_VERSION_MINOR
  integer, parameter :: version_patch = LIBSUPERMESH_VERSION_PATCH
#ifdef LIBSUPERMESH_VERSION_RELEASE
  logical, parameter :: version_release = .true.
#else
  logical, parameter :: version_release = .false.
#endif
  character(len = *), parameter :: version = LIBSUPERMESH_VERSION

end module libsupermesh
