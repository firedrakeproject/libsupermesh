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
module libsupermesh_global_numbering
  ! **********************************************************************
  ! Module to construct the global node numbering map for elements of a
  ! given degree.
  use libsupermesh_adjacency_lists
  use libsupermesh_elements
  use libsupermesh_sparse_tools
  use libsupermesh_fldebug
  use libsupermesh_halo_data_types
  use libsupermesh_halos_allocates
  use libsupermesh_halos_base
!  use halos_debug		! IAKOVOS commented out
  use libsupermesh_halos_numbering
!  use halos_ownership		! IAKOVOS commented out
  use libsupermesh_parallel_tools
  use libsupermesh_linked_lists
  use libsupermesh_mpi_interfaces
  use libsupermesh_fields_base
  
  implicit none

  private
  
! IAKOVOS commented out
!  public :: make_global_numbering_DG, make_boundary_numbering,&
!       & make_global_numbering, element_halo_communicate_visibility, &
!       & make_global_numbering_trace
  public :: make_global_numbering_DG, &
       & make_global_numbering, make_global_numbering_trace

contains

  subroutine make_global_numbering_DG(new_nonods, new_ndglno, Totele,&
       & element, element_halos, new_halos)
    ! Construct a global node numbering for the solution variables in a
    ! Discontinuous Galerkin simulation. This is trivial.
    !
    ! Note that this code is broken for mixed element meshes.
    integer, intent(in) :: totele
    type(element_type), intent(in) :: element

    integer, dimension(:), intent(out) :: new_ndglno
    integer, intent(out) :: new_nonods
    type(halo_type), dimension(:), intent(in), optional :: element_halos
    type(halo_type), dimension(:), intent(out), optional :: new_halos

    integer :: i

! IAKOVOS commented out
    FLAbort("make_global_numbering_DG: Code Commented out")

  end subroutine make_global_numbering_DG
  
  subroutine make_global_numbering_trace(mesh)
    ! Construct a global node numbering for a trace mesh
    !
    ! Note that this code is broken for mixed element meshes.
    type(libsupermesh_mesh_type), intent(inout) :: mesh
    !
    integer :: ele, totele, ni, ele_2, current_global_index
    integer, pointer, dimension(:) :: neigh
    integer :: face_1,face_2,nfaces,i,face_loc, nloc

    totele = mesh%elements
    face_loc = mesh%faces%shape%loc
    nloc = mesh%shape%loc

    !count up how many faces there are
    nfaces = 0
    do ele = 1, totele
       neigh => ele_neigh(mesh,ele)
       do ni = 1, size(neigh)
          ele_2 = neigh(ni)
          if(ele_2<ele) then
             nfaces=nfaces+1
          end if
       end do
    end do
    mesh%nodes = nfaces*face_loc

    !construct mesh%ndglno
    mesh%ndglno = 0
    current_global_index = 0
    do ele = 1, totele
       neigh => ele_neigh(mesh,ele)
       do ni = 1, size(neigh)
          ele_2 = neigh(ni)
          if(ele_2<ele) then
             face_1=ele_face(mesh, ele, ele_2)
             mesh%ndglno((ele-1)*nloc+face_local_nodes(mesh,face_1))&
                  &=current_global_index+(/(i, i=1,face_loc)/)
             if(ele_2>0) then
                !it's not a domain boundary
                !not quite sure how this works in parallel
                face_2=ele_face(mesh, ele_2, ele)
                mesh%ndglno((ele_2-1)*nloc+face_local_nodes(mesh,face_2))&
                     &=current_global_index+(/(i, i=1,face_loc)/)
             end if
             current_global_index = current_global_index + &
                  & mesh%faces%shape%loc
          end if
       end do
    end do
    if(current_global_index /= mesh%nodes) then
       FLAbort('bad global index count in make_global_numbering_trace')
    end if
    if(any(mesh%ndglno==0)) then
       FLAbort('Failed to fully populate trace mesh ndglno')
    end if

  end subroutine make_global_numbering_trace

  subroutine make_global_numbering &
       (new_nonods, new_ndglno, Nonods, Totele, NDGLNO, element, halos,&
       & element_halo, new_halos) 
    ! Construct the global node numbering based on the node numbering for
    ! linear tets given in NDGLNO. 
    integer, intent(in) :: nonods, totele
    integer, intent(in), dimension(:), target :: ndglno
    type(element_type), intent(in) :: element
    !! The level 1 and 2 halos associated with the incoming mesh.
    type(halo_type), dimension(:), intent(in), optional :: halos
    !! The full element halo associated with these meshes.
    type(halo_type), intent(in), optional :: element_halo
    !! The level 1 and 2 halos associated with the new mesh.
    type(halo_type), dimension(:), intent(out), optional :: new_halos

    integer, dimension(:), intent(out) :: new_ndglno
    integer, intent(out) :: new_nonods

    ! Adjacency lists.
    type(csr_sparsity) :: NEList, NNList, EEList

    logical :: D3
    integer :: dim, nloc, snloc, faces, owned_nodes

    ! Number of nodes associated with each object.
    integer :: face_len, edge_len, element_len

    ! Total nodes associated with an element
    integer :: element_tot_len

    integer :: ele, ele2, node, node2, new_node, new_node2, i, j, k, halo

    integer, dimension(:), allocatable :: n

    ! Process number of this processor
    integer :: rank

    ! Owner and halo_level of all the nodes
    integer, dimension(:), allocatable :: node_owner, receive_halo_level, &
         & new_receive_halo_level
    integer, dimension(:,:), allocatable :: new_node_owner
    integer :: this_node_owner, this_receive_halo_level, this_node_winner

    ! Cache for positions in ndglno which are currently being worked on and
    ! new values to go in them
    integer, dimension(:), allocatable :: ndglno_pos, ndglno_val, face_nodes

    ! Map from old node numbers to new ones.
    integer, dimension(:), allocatable :: node_map
    
    ! In each case, halos is the last dimension of halos.
    type(ilist), dimension(:), allocatable :: this_send_targets, node_targets,&
         & node2_targets
    type(ilist), dimension(:,:), allocatable :: old_send_targets, new_send_targets
    type(inode), pointer :: this_target

    ! Nodes in current element
    integer, pointer, dimension(:) :: ele_node

    ! Flag for whether halos are being calculated:
    logical :: have_halos

! IAKOVOS commented out
    FLAbort("make_global_numbering: Code Commented out")

  end subroutine make_global_numbering

end module libsupermesh_global_numbering
