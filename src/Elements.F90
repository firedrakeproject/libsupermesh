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
module libsupermesh_elements
  !!< This module provides derived types for finite elements and associated functions.
  use libsupermesh_element_numbering
  use libsupermesh_quadrature
  use libsupermesh_FLDebug
  use libsupermesh_polynomials
  use libsupermesh_reference_counting
  implicit none

  type element_type
     !!< Type to encode shape and quadrature information for an element.
     integer :: dim !! 2d or 3d?
     integer :: loc !! Number of nodes.
     integer :: ngi !! Number of gauss points.
     integer :: degree !! Polynomial degree of element.
     !! Shape functions: n is for the primitive function, dn is for partial derivatives, dn_s is for partial derivatives on surfaces. 
     !! n is loc x ngi, dn is loc x ngi x dim
     !! dn_s is loc x ngi x face x dim 
     real, pointer :: n(:,:)=>null(), dn(:,:,:)=>null()
     real, pointer :: n_s(:,:,:)=>null(), dn_s(:,:,:,:)=>null()
     !! Polynomials defining shape functions and their derivatives.
     type(polynomial), dimension(:,:), pointer :: spoly=>null(), dspoly=>null()
     !! Link back to the node numbering used for this element.
     type(ele_numbering_type), pointer :: numbering=>null()
     !! Link back to the quadrature used for this element.
     type(quadrature_type) :: quadrature
     type(quadrature_type), pointer :: surface_quadrature=>null()
     !! Pointer to the superconvergence data for this element.
     type(superconvergence_type), pointer :: superconvergence=>null()
     !! Pointer to constraints data for this element
     type(constraints_type), pointer :: constraints=>null()
     !! Reference count to prevent memory leaks.
     type(refcount_type), pointer :: refcount=>null()
     !! Dummy name to satisfy reference counting
     character(len=0) :: name
  end type element_type
  
  type superconvergence_type
    !!< A structure to represent the superconvergent points of the element in question.
    !!< This is in this module because it has to be in element_type,
    !!< but Superconvergence.F90 depends on Elements.F90. So Elements.F90
    !!< cannot depend on Superconvergence.F90. (Fortran is a real pain.)
    !! Number of superconvergent points
    integer :: nsp 
    !! Locations of superconvergent points in local coordinates
    !! allocated to nsp x loc
    real, pointer :: l(:, :)
    !! Shape functions at each superconvergent point.
    !! loc x nsp
    real, pointer :: n(:, :)
    !! Derivatives of shape functions at each superconvergent point
    !! loc x nsp x ndim
    real, pointer :: dn(:, :, :)
  end type superconvergence_type
  
  type constraints_type
     !!< A type to encode the constraints from the local Lagrange basis for 
     !!< (Pn)^d vector-valued elements to another local basis, possibly for 
     !!< a proper subspace. This new basis must have DOFs consisting
     !!< of either normal components on faces corresponding to a Lagrange
     !!< basis for the normal component when restricted to each face,
     !!< or coefficients of basis
     !!< functions with vanishing normal components on all faces.
     !! type of constraints
     integer :: type
     !! local dimension
     integer :: dim
     !! order of local Lagrange basis
     integer :: degree
     !! number of nodes for local Lagrange basis
     integer :: loc
     !! Number of constraints
     integer :: n_constraints
     !! basis of functions that are orthogonal to the 
     !! constrained vector space 
     !! dimension n_constraints x loc x dim
     real, pointer :: orthogonal(:,:,:)=> null()
  end type constraints_type
  
  integer, parameter :: CONSTRAINT_NONE =0, CONSTRAINT_BDFM = 1,&
       & CONSTRAINT_RT = 2, CONSTRAINT_BDM = 3
       
  interface deallocate
     module procedure deallocate_element
     module procedure deallocate_constraints
  end interface

  interface local_coords
     module procedure element_local_coords
  end interface

  interface eval_shape
    module procedure eval_shape_node, eval_shape_all_nodes
  end interface
  
  interface eval_dshape
    module procedure eval_dshape_node, eval_dshape_all_nodes
  end interface
 
#include "Reference_count_interface_element_type.F90"
  
contains

  subroutine deallocate_element(element, stat)
    type(element_type), intent(inout) :: element
    integer, intent(out), optional :: stat
    
    integer :: lstat, tstat
    integer :: i,j

    tstat = 0
    lstat = 0

    call decref(element)
    if (has_references(element)) then
       ! There are still references to this element so we don't deallocate.
       return
    end if

    call deallocate(element%quadrature)

    if(associated(element%spoly)) then
      do i=1,size(element%spoly,1)
        do j=1,size(element%spoly,2)
            call deallocate(element%spoly(i,j))
        end do
      end do
      deallocate(element%spoly, stat=tstat)
    end if
    lstat=max(lstat,tstat)

    if (associated(element%n_s)) deallocate(element%n_s,element%dn_s)

    if(associated(element%dspoly)) then
      do i=1,size(element%dspoly,1)
        do j=1,size(element%dspoly,2)
            call deallocate(element%dspoly(i,j))
        end do
      end do
      deallocate(element%dspoly, stat=tstat)
    end if
    lstat=max(lstat,tstat)

    deallocate(element%n,element%dn, stat=tstat)
    lstat=max(lstat,tstat)

    if(associated(element%constraints)) then
       call deallocate(element%constraints,stat=tstat)
       lstat = max(lstat,tstat)

       deallocate(element%constraints, stat=tstat)
       lstat = max(lstat,tstat)
    end if
    if (present(stat)) then
       stat=lstat
    else if (lstat/=0) then
!       FLAbort("Unable to deallocate element.")		! ToDo
    end if

  end subroutine deallocate_element
  
  subroutine deallocate_constraints(constraint, stat)
    type(constraints_type), intent(inout) :: constraint
    integer, intent(out), optional :: stat
    
    integer :: lstat

    lstat = 0

    if(associated(constraint%orthogonal)) then
       deallocate(constraint%orthogonal,stat=lstat)
    end if

    if (present(stat)) then
       stat=lstat
    else if (lstat/=0) then
!       FLAbort("Unable to deallocate constraints.")		! ToDo
    end if

  end subroutine deallocate_constraints
  
  function element_local_coords(n, element) result (coords)
    !!< Work out the local coordinates of node n in element. This is just a
    !!< wrapper function which allows local_coords to be called on an element
    !!< instead of on an element numbering.
    integer, intent(in) :: n
    type(element_type), intent(in) :: element    
    real, dimension(size(element%numbering%number2count, 1)) :: coords
    
    coords=local_coords(n, element%numbering)

  end function element_local_coords

  pure function eval_shape_node(shape, node,  l) result(eval_shape)
    ! Evaluate the shape function for node node local coordinates l
    real :: eval_shape
    type(element_type), intent(in) :: shape
    integer, intent(in) :: node
    real, dimension(size(shape%spoly,1)), intent(in) :: l

    integer :: i

    eval_shape=1.0
          
    do i=1,size(shape%spoly,1)
       
       ! Raw shape function
       eval_shape=eval_shape*eval(shape%spoly(i,node), l(i))
             
    end do

  end function eval_shape_node

  pure function eval_shape_all_nodes(shape, l) result(eval_shape)
    ! Evaluate the shape function for all locations at local coordinates l
    type(element_type), intent(in) :: shape
    real, dimension(size(shape%spoly,1)), intent(in) :: l
    real, dimension(shape%loc) :: eval_shape

    integer :: i,j

    eval_shape=1.0

    do j=1,shape%loc

      do i=1,size(shape%spoly,1)

        ! Raw shape function
        eval_shape(j)=eval_shape(j)*eval(shape%spoly(i,j), l(i))

      end do

    end do

  end function eval_shape_all_nodes
  
  pure function eval_dshape_node(shape, node,  l) result(eval_dshape)
    !!< Evaluate the derivatives of the shape function for location node at local
    !!< coordinates l 
    type(element_type), intent(in) :: shape
    integer, intent(in) :: node
    real, dimension(:), intent(in) :: l
    real, dimension(shape%dim) :: eval_dshape

    select case(shape%numbering%family)
       
    case (FAMILY_SIMPLEX)

       eval_dshape=eval_dshape_simplex(shape, node,  l)

    case (FAMILY_CUBE)

       eval_dshape=eval_dshape_cube(shape, node,  l)

    case default
       ! Invalid element family. Return a really big number to stuff things
       ! quickly. 

       eval_dshape=huge(0.0)

    end select
    
  end function eval_dshape_node

  function eval_dshape_all_nodes(shape, l) result(eval_dshape)
    type(element_type), intent(in) :: shape
    real, dimension(:), intent(in) :: l
    real, dimension(shape%loc, shape%dim) :: eval_dshape

    integer :: loc

    do loc=1,shape%loc
      eval_dshape(loc, :) = eval_dshape_node(shape, loc, l)
    end do
  end function eval_dshape_all_nodes
  
  pure function eval_dshape_simplex(shape, loc,  l) result (eval_dshape)
    !!< Evaluate the derivatives of the shape function for location loc at local
    !!< coordinates l 
    !!<
    !!< This version of the function applies to members of the simplex
    !!< family including the interval.
    type(element_type), intent(in) :: shape
    integer, intent(in) :: loc
    real, dimension(shape%dim+1), intent(in) :: l
    real, dimension(shape%dim) :: eval_dshape
    
    integer :: i,j
    ! Derivative of the dependent coordinate with respect to the other
    ! coordinates:
    real, dimension(shape%dim) :: dl4dl

    ! Find derivative of dependent coordinate
    dl4dl=diffl4(shape%numbering%vertices, shape%dim)

    do i=1,shape%dim
       ! Directional derivatives.
       
       ! The derivative has to take into account the dependent
       ! coordinate. In 3D:
       !
       !  S=P1(L1)P2(L2)P3(L3)P4(L4)
       !
       !  dS        / dP1     dL4 dP4  \
       !  --- = P2P3| ---P4 + ---*---P1|
       !  dL1       \ dL1     dL1 dL4  /
       !
       
       ! Expression in brackets.
       eval_dshape(i)=eval(shape%dspoly(i,loc), l(i))&
            *eval(shape%spoly(shape%dim+1,loc),l(shape%dim+1))&
            + dl4dl(i)&
            *eval(shape%dspoly(shape%dim+1,loc), l(shape%dim+1)) &
            *eval(shape%spoly(i,loc),l(i))
             
       ! The other terms
       do j=1,shape%dim
          if (j==i) cycle
          
          eval_dshape(i)=eval_dshape(i)*eval(shape%spoly(j,loc), l(j))
       end do
       
    end do

  end function eval_dshape_simplex
  
  pure function eval_dshape_cube(shape, loc,  l) result (eval_dshape)
    !!< Evaluate the derivatives of the shape function for location loc at local
    !!< coordinates l 
    !!<
    !!< This version of the function applies to members of the hypercube
    !!< family. Note that this does NOT include the interval.
    type(element_type), intent(in) :: shape
    integer, intent(in) :: loc
    real, dimension(shape%dim+1), intent(in) :: l
    real, dimension(shape%dim) :: eval_dshape

    integer :: i,j

    do i=1,shape%dim
       eval_dshape(i)=1.0
       ! Directional derivatives.
       do j=1,shape%dim
          if(i==j) then
            eval_dshape(i)=eval_dshape(i)*eval(shape%dspoly(j,loc), l(j))
          else
            eval_dshape(i)=eval_dshape(i)*eval(shape%spoly(j,loc), l(j))
          end if          
       end do
    
    end do

  end function eval_dshape_cube
  
  pure function diffl4(vertices, dimension)
    ! Derivative of the dependent coordinate with respect to the other
    ! coordinates. 
    integer, intent(in) :: vertices, dimension
    real, dimension(dimension) :: diffl4

    if (vertices==dimension+1) then
       ! Simplex. Dependent coordinate depends on all other coordinates. 
       diffl4=-1.0
       
    else if (vertices==2**dimension) then
       ! Hypercube. The dependent coordinate is redundant.
       diffl4=0.0
    
    else if (vertices==6.and.dimension==3) then
       ! Wedge. First coordinate is independent.
       diffl4=(/0.0,-1.0,-1.0/)

    else
       ! No output permitted in a pure procedure so we return a big number to stuff
       ! things up quickly.
       diffl4=huge(0.0)
    end if
       
  end function diffl4
  
#include "Reference_count_element_type.F90"

end module libsupermesh_elements
