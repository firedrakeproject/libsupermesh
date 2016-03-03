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
  integer, parameter :: version_release = libsupermesh_VERSION_RELEASE
  character(len = *), parameter :: version = libsupermesh_VERSION_STRING

end module libsupermesh
