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
module libsupermesh_transform_elements
  ! Module to calculate element transformations from local to physical
  ! coordinates.
  use libsupermesh_quadrature
  use libsupermesh_elements
  use libsupermesh_vector_tools
  use libsupermesh_parallel_tools, only: abort_if_in_parallel_region
  use libsupermesh_fields_base
!  use cv_faces, only: cv_faces_type	! IAKOVOS commented out
  use libsupermesh_eventcounter
!  use memory_diagnostics		! IAKOVOS commented out
  
  implicit none
  
  interface transform_facet_to_physical
    module procedure transform_facet_to_physical_full
  end interface transform_facet_to_physical
  
  interface retrieve_cached_face_transform
     module procedure retrieve_cached_face_transform_full
  end interface
  
  private
  public :: transform_facet_to_physical
  
  integer, parameter :: cyc3(1:5)=(/ 1, 2, 3, 1, 2 /)  
  
  logical, save :: cache_transform_elements=.true.
  
  real, dimension(:,:), allocatable, save :: face_normal_cache
  real, dimension(:), allocatable, save :: face_detJ_cache
  ! Record which element is on the other side of the last n/2 elements.
  integer, dimension(:), allocatable, save :: face_cache
  
  ! The reference count id of the positions mesh being cached.
  integer, save :: face_position_id=-1
  integer, save :: face_last_mesh_movement=-1
  
contains
  function retrieve_cached_face_transform_full(X, face, &
       & normal_local, detJ_local) result (cache_valid)
    !!< Determine whether the transform cache is valid for this operation.
    !!< 
    !!< If caching is applicable and the cache is not ready, set up the
    !!< cache and then return true.
    type(vector_field), intent(in) :: X
    integer, intent(in) :: face
    !! Face determinant
    real, intent(out) :: detJ_local
    !! Face normal
    real, dimension(X%dim), intent(out) :: normal_local

    logical :: cache_valid

    cache_valid=.true.
    
    if (X%refcount%id/=face_position_id) then
!       ewrite(2,*) "Reference count identity of X has changed."
       cache_valid=.false.
       if (X%name/="Coordinate") then
          !!< If Someone is not calling this on the main Coordinate field
          !!< then we're screwed anyway.
          return
       end if
    
    else if(eventcount(EVENT_MESH_MOVEMENT)/=face_last_mesh_movement) then
!       ewrite(2,*) "Mesh has moved."
       cache_valid=.false.

    end if
       
    if (.not.cache_valid) then
       call construct_face_cache(X)
       cache_valid=.true.
    end if

    detJ_local=face_detJ_cache(abs(face_cache(face)))  
    normal_local=sign(1,face_cache(face))*face_normal_cache(:,abs(face_cache(face)))
    
  end function retrieve_cached_face_transform_full

  subroutine construct_face_cache(X)
    !!< The cache is invalid so make a new one.
    type(vector_field), intent(in) :: X
    
    integer :: elements, ele, i, current_face, face, face2, faces, n,&
         & unique_faces
    !! Note that if X is not all linear simplices we are screwed. 
    real, dimension(X%dim, ele_loc(X,1)) :: X_val
    real, dimension(X%dim, face_loc(X,1)) :: X_f
    real, dimension(X%dim, X%dim-1) :: J
    type(element_type), pointer :: X_shape_f
    real :: detJ
    integer, dimension(:), pointer :: neigh
    
!    ewrite(1,*) "Reconstructing element geometry cache."

    call abort_if_in_parallel_region

    face_position_id=X%refcount%id
    face_last_mesh_movement=eventcount(EVENT_MESH_MOVEMENT)

    if (allocated(face_detJ_cache)) then
#ifdef HAVE_MEMORY_STATS
       call register_deallocation("transform_cache", "real", &
            & size(face_detJ_cache)+size(face_normal_cache))
       call register_deallocation("transform_cache", "integer", &
            & size(face_cache))
#endif
       deallocate(face_detJ_cache, face_normal_cache, face_cache)
    end if

    elements=element_count(X)
    faces=face_count(X)
    !! This counts 1/2 for each interior face and 1 for each surface face.
    unique_faces=unique_face_count(X%mesh)

    allocate(face_detJ_cache(unique_faces), &
         face_normal_cache(X%dim,unique_faces), &
         face_cache(faces))
#ifdef HAVE_MEMORY_STATS
    call register_allocation("transform_cache", "real", &
         & size(face_detJ_cache)+size(face_normal_cache))
    call register_allocation("transform_cache", "integer", &
         & size(face_cache))
#endif

    current_face=0
    do ele=1,elements
       neigh=>ele_neigh(X, ele)
       X_val=ele_val(X,ele)

       do n=1,size(neigh)

          if (neigh(n)<0) then
             face=ele_face(X, ele, neigh(n))

             current_face=current_face+1
             face_cache(face)=current_face
             
          else
                       
             face=ele_face(X, ele, neigh(n))
             face2=ele_face(X, neigh(n), ele)
             
             ! Only do this once for each face pairl
             if (face>face2) then
                cycle
             end if

             current_face=current_face+1
             
             face_cache(face)=current_face
             face_cache(face2)=-current_face
             
          end if
          
          
          X_f=face_val(X,face)
          X_shape_f=>face_shape(X,face)

          !     |- dx  dx  -|
          !     |  dL1 dL2  |
          !     |           |
          !     |  dy  dy   |
          ! J = |  dL1 dL2  |
          !     |           |
          !     |  dz  dz   |
          !     |- dL1 dL2 -|
          
          ! Form Jacobian.
          J=matmul(X_f(:,:), x_shape_f%dn(:, 1, :))

          detJ=0.0
          ! Calculate determinant.
          select case (X%dim)
          case(1)
             detJ=1.0
          case(2)
             detJ = sqrt(J(1,1)**2 + J(2,1)**2)
          case(3)
             do i=1,3
                detJ=detJ+ &
                     (J(cyc3(i+2),1)*J(cyc3(i+1),2)-J(cyc3(i+2),2)*J(cyc3(i+1),1))**2
             end do
             detJ=sqrt(detJ)
          case default
             FLAbort("Unsupported dimension specified.  Universe is 3 dimensional (sorry Albert).")   
          end select

          ! Calculate normal.
          face_normal_cache(:,current_face)=normgi(X_val,X_f,J)
          face_detJ_cache(current_face)=detJ
   
       end do
    end do
    assert(current_face==unique_faces)

  end subroutine construct_face_cache
  
  function unique_face_count(mesh) result (face_count)
    !!< Count the number of geometrically unique faces in mesh,
    type(mesh_type), intent(in) :: mesh
    integer :: face_count

    integer :: ele
    integer, dimension(:), pointer :: neigh

    face_count=0

    do ele=1,element_count(mesh)
       neigh=>ele_neigh(mesh, ele)
       
       ! Count 1 for interior and 2 for surface.
       face_count=face_count+sum(merge(1,2,neigh>0)) 

    end do

    face_count=(face_count+1)/2
    
  end function unique_face_count

  subroutine transform_facet_to_physical_full(X, face, detwei_f, normal)

    ! Coordinate transformations for facet integrals. 
    ! Calculate the transformed quadrature
    ! weights as a side bonus.
    !
    ! For facet integrals, we also need to know the facet outward
    ! pointing normal.
    !
    ! In this case it is only the determinant of the Jacobian which is
    ! required.

    ! Column n of X is the position of the nth node of the adjacent element
    ! (this is only used to work out the orientation of the boundary)
    type(vector_field), intent(in) :: X
    ! The face to transform.
    integer, intent(in) :: face
    ! Quadrature weights for physical coordinates for integration over the boundary.
    real, dimension(:), intent(out), optional :: detwei_f
    ! Outward normal vector. (dim x x_shape_f%ngi)
    real, dimension(:,:), intent(out) :: normal

    
    ! Column n of X_f is the position of the nth node on the facet.
    real, dimension(X%dim,face_loc(X,face)) :: X_f
    ! Column n of X_f is the position of the nth node on the facet.
    real, dimension(X%dim,ele_loc(X,face_ele(X,face))) :: X_val
    ! shape function coordinate interpolation on the boundary
    type(element_type), pointer :: x_shape_f


    ! Jacobian matrix and its inverse.
    real, dimension(X%dim,mesh_dim(X)-1) :: J
    ! Determinant of J
    real :: detJ
    ! Whether the cache can be used
    logical :: cache_valid
    

    integer :: gi, i, compute_ngi

    x_shape_f=>face_shape(X,face)

#ifdef DDEBUG
    assert(size(normal,1)==X%dim)
#endif
#ifdef DDEBUG
    if (present(detwei_f)) then
       assert(size(detwei_f)==x_shape_f%ngi)
    end if
#endif
    
    if (.not.(x_shape_f%degree==1 .and. x_shape_f%numbering%family==FAMILY_SIMPLEX)) then
      ! for non-linear compute on all gauss points
      compute_ngi=x_shape_f%ngi
      cache_valid=.false.
    else
      ! for linear: compute only the first and copy the rest
      if (cache_transform_elements) then
         cache_valid=retrieve_cached_face_transform(X, face, normal(:,1),&
              & detJ)
      else
         cache_valid=.false.
      end if
      if (cache_valid) then
         compute_ngi=0
      else
         compute_ngi=1
      end if
      
    end if

    if (.not.cache_valid) then
       X_val=ele_val(X, face_ele(X,face))
       X_f=face_val(X, face)
    end if

    ! Loop over quadrature points.
    quad_loop: do gi=1, compute_ngi

       !     |- dx  dx  -|
       !     |  dL1 dL2  |
       !     |           |
       !     |  dy  dy   |
       ! J = |  dL1 dL2  |
       !     |           |
       !     |  dz  dz   |
       !     |- dL1 dL2 -|

       ! Form Jacobian.
       J=matmul(X_f(:,:), x_shape_f%dn(:, gi, :))

       detJ=0.0
       ! Calculate determinant.
       select case (mesh_dim(X))
       case(1)
          detJ=1.0
       case(2)
          select case (X%dim)
          case(2)
             detJ = sqrt(J(1,1)**2 + J(2,1)**2)
          case(3)
             detJ = sqrt(sum(J(:,1)**2))
          case default
             FLAbort("Unsupported dimension specified")
          end select
       case(3)
          select case (X%dim)
          case(3)
             do i=1,3
                detJ=detJ+ &
                     (J(cyc3(i+2),1)*J(cyc3(i+1),2)-J(cyc3(i+2),2)*J(cyc3(i&
                     &+1),1))**2
             end do
             detJ=sqrt(detJ)
          case default
             FLAbort("Unsupported dimension specified")
          end select
       case default
          FLAbort("Unsupported dimension specified.  Universe is 3 dimensional (sorry Albert).")   
       end select

       ! Calculate transformed quadrature weights.
       if(present(detwei_f)) then
          detwei_f(gi)=detJ*x_shape_f%quadrature%weight(gi)
       end if
       ! Calculate normal.
       normal(:,gi)=normgi(X_val,X_f,J)

    end do quad_loop
      
    ! copy the value at gi==1 to the rest of the gauss points
    if(present(detwei_f)) then
      do gi=compute_ngi+1, x_shape_f%ngi
         ! uses detJ from above
         detwei_f=detJ*x_shape_f%quadrature%weight
      end do
    end if

    do gi=compute_ngi+1, x_shape_f%ngi
       normal(:,gi)=normal(:,1)
    end do    
    
  end subroutine transform_facet_to_physical_full
  
  function NORMGI(X, X_f, J)
    ! Calculate the normal at a given quadrature point,
    real, dimension(:,:), intent(in) :: J
    real, dimension(size(J,1)) :: normgi
    ! Element and normal node locations respectively.
    real, dimension (:,:), intent(in) :: X, X_f
    ! Facet Jacobian.

    ! Outward pointing not necessarily normal vector.
    real, dimension(3) :: outv

    integer :: ldim

    ldim = size(J,1)

    ! Outv is the vector from the element centroid to the facet centroid.
    outv(1:ldim) = sum(X_f,2)/size(X_f,2)-sum(X,2)/size(X,2)

    select case (ldim)
    case(1)
       normgi = 1.0
    case (2)
       normgi = (/ -J(2,1), J(1,1) /)
    case (3)
       normgi=cross_product(J(:,1),J(:,2))
    case default
       FLAbort("Unsupported dimension specified.  Universe is 3 dimensional (sorry Albert).")   
    end select

    ! Set correct orientation.
    normgi=normgi*dot_product(normgi, outv(1:ldim) )

    ! normalise
    normgi=normgi/sqrt(sum(normgi**2))

  contains

    function cross_product(vector1,vector2) result (prod)
      real, dimension(3) :: prod
      real, dimension(3), intent(in) :: vector1, vector2

      integer :: i     

      forall(i=1:3)
         prod(i)=vector1(cyc3(i+1))*vector2(cyc3(i+2))&
              -vector1(cyc3(i+2))*vector2(cyc3(i+1))
      end forall

    end function cross_product

  end function NORMGI
  
end module libsupermesh_transform_elements
