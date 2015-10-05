#include "fdebug.h"

module libsupermesh_mesh_stitching
  
  use libsupermesh_fldebug
  
  implicit none

#include "mpif.h"
  
  private
  
  public :: stich_meshes
  
contains  

  subroutine stich_meshes(positions, enlist, un, positions_add, enlist_add, un_add, positions_new, enlist_new)
    ! dim x nnodes
    real, dimension(:, :), intent(in) :: positions
    ! loc x nelements
    integer, dimension(:, :), intent(in) :: enlist
    ! nnodes
    integer, dimension(:) :: un
    ! dim x nnodes_add
    real, dimension(:, :), intent(in) :: positions_add
    ! loc x nelements_add
    integer, dimension(:, :), intent(in) :: enlist_add
    ! nnodes_add
    integer, dimension(:) :: un_add
    ! dim x nnodes_new
    real, dimension(:, :), allocatable, intent(out) :: positions_new
    ! loc x nelements_new
    integer, dimension(:, :), allocatable, intent(out) :: enlist_new
  
  end subroutine stich_meshes    

end module libsupermesh_mesh_stitching
