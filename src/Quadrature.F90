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
module libsupermesh_quadrature
  !!< This module implements quadrature of varying degrees for a number of
  !!< elements. Quadrature information is used to numerically evaluate
  !!< integrals over an element. 
  use libsupermesh_FLDebug
  use libsupermesh_reference_counting
!  use wandzura_quadrature		! IAKOVOS commented out
!  use grundmann_moeller_quadrature	! IAKOVOS commented out
  use libsupermesh_vector_tools
  implicit none

  private
  
  type generator_type
     !!< The generator type is an encoding of a quadrature generator of the
     !!< type used in the encyclopedia of cubature. This type is only used
     !!< internally in the quadrature module.
     integer, dimension(:,:), pointer :: permutation
     real, dimension(:), pointer :: coords
     real :: weight
  end type generator_type

  type quadrature_template
     !!< A data type which defines a quadrature rule. These are only
     !!< directly used inside the quadrature module.

     !! A quadrature is defined by a set of generators.
     type(generator_type), dimension(:), pointer :: generator
     !! Dimension of the space we are in and the degree of accuracy of the
     !! quadrature. 
     integer :: dim, degree
     !! Ngi is number of quadrature points. Vertices is number of vertices. These
     !! names are chosen for consistency with the rest of fluidity.
     integer :: ngi, vertices
  end type quadrature_template
  
  type quadrature_type
     !!< A data type which describes quadrature information. For most
     !!< developers, quadrature can be treated as an opaque data type which
     !!< will only be encountered when creating element_type variables to
     !!< represent shape functions.  
     integer :: dim !! Dimension of the elements for which quadrature
     !!< is required.  
     integer :: degree !! Degree of accuracy of quadrature. 
     integer :: vertices !! Number of vertices of the element.
     integer :: ngi !! Number of quadrature points.
     real, pointer :: weight(:)=>null() !! Quadrature weights.
     real, pointer :: l(:,:)=>null() !! Locations of quadrature points.
     character(len=0) :: name !! Fake name for reference counting.
     !! Reference count to prevent memory leaks.
     type(refcount_type), pointer :: refcount=>null()
     integer :: family
  end type quadrature_type
  
  type(quadrature_template), dimension(8), target, save, public :: tet_quads
  type(quadrature_template), dimension(8), target, save, public :: tri_quads
  type(quadrature_template), dimension(8), target, save, public :: interval_quads
  type(quadrature_template), dimension(6), target, save, public :: hex_quads
  type(quadrature_template), dimension(6), target, save, public :: quad_quads
  type(quadrature_template), dimension(1), target, save, public :: point_quad
  
  character(len=100), save, public :: quadrature_error_message=""
  
  !! Unsupported vertex count.
  integer, parameter, public :: QUADRATURE_VERTEX_ERROR=1 
  !! Quadrature degree requested is not available.
  integer, parameter, public :: QUADRATURE_DEGREE_ERROR=2 
  !! Elements with this number of dimensions are not available.
  integer, parameter, public :: QUADRATURE_DIMENSION_ERROR=3 
  !! Unsupported number of quadrature points.
  integer, parameter, public :: QUADRATURE_NGI_ERROR=4 
  !! Not enough arguments specified.
  integer, parameter, public :: QUADRATURE_ARGUMENT_ERROR=5
  
  integer, parameter :: FAMILY_COOLS=0, FAMILY_WANDZURA=1, FAMILY_GM=2

  interface allocate
     module procedure allocate_quad
  end interface
  
  interface deallocate
     module procedure deallocate_quad
  end interface
  
#include "Reference_count_interface_quadrature_type.F90"
  
  public make_quadrature, deallocate, quadrature_type, &
       & operator(==), incref, addref, decref

contains

  function make_quadrature(vertices, dim, degree, ngi, family, stat) result (quad)
    !!< Given information about a quadrature, return a quad type encoding
    !!< that quadrature.
    type(quadrature_type) :: quad
    !! Using vertices and dimension it is possible to determine what shape we are
    !! using. At this stage we assume that no-one will require elements in
    !! the shape of the tetragonal antiwedge!
    integer, intent(in) :: vertices, dim
    !! Ngi is the old way of specifying quadrature. This should really be
    !! done via degree. At least one of these must be specified. If both are
    !! specified then ngi is used.
    integer, intent(in), optional :: degree, ngi
    !! Which family of quadrature you'd like to use.
    integer, intent(in), optional :: family
    !! Status argument - zero for success non-zero otherwise.
    integer, intent(out), optional :: stat

    ! The set of quadrature templates for this shape of element.
    type(quadrature_template), dimension(:), pointer :: template_set
    ! The quadrature template we will use.
    type(quadrature_template), pointer :: template
    ! Number of local coordinates
    integer coords

    integer :: lfamily
    integer :: wandzura_rule_idx, wandzura_rule_degree, max_wandzura_rule, wandzura_order
    real, dimension(2, 3) :: wandzura_ref_tri
    real, dimension(3, 3) :: wandzura_ref_map
    real, dimension(:, :), allocatable :: tmp_coordinates
    integer :: gi

    integer :: gm_rule, gm_order, vertex
    real, dimension(:, :), allocatable :: gm_ref_simplex
    real, dimension(:, :), allocatable :: gm_ref_map

! IAKOVOS commented out
    FLAbort("make_quadrature: Code Commented out")

  end function make_quadrature
  
  subroutine allocate_quad(quad, vertices, ngi, coords, stat)
    !!< Allocate memory for a quadrature type. Note that this is done
    !!< automatically in make_quadrature. 
    type(quadrature_type), intent(inout) :: quad
    !! Vertices is the number of vertices. Ngi is the number of quadrature
    !! points. Coords the number of local coords
    integer, intent(in) :: vertices, ngi, coords
    !! Stat returns zero for successful completion and nonzero otherwise.
    integer, intent(out), optional :: stat
  
    integer :: lstat

    allocate(quad%weight(ngi), quad%l(ngi,coords), stat=lstat)
    
    quad%vertices=vertices
    quad%ngi=ngi

    nullify(quad%refcount) ! Hack for gfortran component initialisation
    !                         bug.

    call addref(quad)

    if (present(stat)) then
       stat=lstat
    else if (lstat/=0) then
       FLAbort("Error allocating quad")
    end if

  end subroutine allocate_quad

  subroutine deallocate_quad(quad,stat)
    !!< Since quadrature types contain pointers it is necessary to
    !!< explicitly deallocate them.
    !! The quadrature type to be deallocated.
    type(quadrature_type), intent(inout) :: quad
    !! Stat returns zero for successful completion and nonzero otherwise.
    integer, intent(out), optional :: stat

    integer :: lstat

    call decref(quad)
    if (has_references(quad)) then
       ! There are still references to this quad so we don't deallocate.
       return
    end if
    
    deallocate(quad%weight,quad%l, stat=lstat)
    
    if (present(stat)) then
       stat=lstat
    else if (lstat/=0) then
       FLAbort("Error deallocating quad")
    end if

  end subroutine deallocate_quad
  
#include "Reference_count_quadrature_type.F90"

  recursive function factorial(n) result(f)
    ! Calculate n!
    integer :: f
    integer, intent(in) :: n

    if (n==0) then
       f=1
    else
       f=n*factorial(n-1)
    end if

  end function factorial
  
end module libsupermesh_quadrature
