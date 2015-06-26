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
  
  interface transform_to_physical
!    module procedure transform_to_physical_full, transform_to_physical_detwei
    module procedure transform_to_physical_detwei
  end interface
  
  interface transform_facet_to_physical
    module procedure transform_facet_to_physical_full
  end interface transform_facet_to_physical
  
  interface retrieve_cached_transform
!     module procedure retrieve_cached_transform_full, &
!          retrieve_cached_transform_det
     module procedure retrieve_cached_transform_det
  end interface

  interface retrieve_cached_face_transform
     module procedure retrieve_cached_face_transform_full
  end interface
  
  private
  public :: transform_to_physical, transform_facet_to_physical
  
  integer, parameter :: cyc3(1:5)=(/ 1, 2, 3, 1, 2 /)  
  
  logical, save :: cache_transform_elements=.true.
  real, dimension(:,:,:), allocatable, save :: invJ_cache
  real, dimension(:,:,:), allocatable, save :: J_T_cache
  real, dimension(:), allocatable, save :: detJ_cache

  real, dimension(:,:), allocatable, save :: face_normal_cache
  real, dimension(:), allocatable, save :: face_detJ_cache
  ! Record which element is on the other side of the last n/2 elements.
  integer, dimension(:), allocatable, save :: face_cache
  
  ! The reference count id of the positions mesh being cached.
  integer, save :: position_id=-1
  integer, save :: last_mesh_movement=-1
  integer, save :: face_position_id=-1
  integer, save :: face_last_mesh_movement=-1
  
contains

  function retrieve_cached_transform_det(X, ele, detJ_local) &
       result (cache_valid)
    !!< Determine whether the transform cache is valid for this operation.
    !!< 
    !!< If caching is applicable and the cache is not ready, set up the
    !!< cache and then return true.
    type(vector_field), intent(in) :: X
    integer, intent(in) :: ele
    !! Local version of the determinant of J
    real, intent(out) :: detJ_local

    logical :: cache_valid
    
    cache_valid=.true.

    if (X%refcount%id/=position_id) then
       cache_valid=.false.
       if (X%name/="Coordinate") then
          !!< If Someone is not calling this on the main Coordinate field
          !!< then we're screwed anyway.
          return
       end if
       
!       ewrite(2,*) "Reference count identity of X has changed."
    else if(eventcount(EVENT_MESH_MOVEMENT)/=last_mesh_movement) then
!       ewrite(2,*) "Mesh has moved."
       cache_valid=.false.

    end if
       
    if (.not.cache_valid) then
       call construct_cache(X)
       cache_valid=.true.
    end if
    
    detJ_local=detJ_cache(ele)    
    
  end function retrieve_cached_transform_det
  
  subroutine construct_cache(X)
    !!< The cache is invalid so make a new one.
    type(vector_field), intent(in) :: X
    
    integer :: elements, ele, i, k
    !! Note that if X is not all linear simplices we are screwed. 
    real, dimension(X%dim, ele_loc(X,1)) :: X_val
    type(element_type), pointer :: X_shape

!    ewrite(1,*) "Reconstructing element geometry cache."

    call abort_if_in_parallel_region

    position_id=X%refcount%id
    last_mesh_movement=eventcount(EVENT_MESH_MOVEMENT)

    if (allocated(invJ_cache)) then
#ifdef HAVE_MEMORY_STATS
       call register_deallocation("transform_cache", &
            "real", size(invJ_cache)+size(J_T_cache)+size(detJ_cache))
#endif
       deallocate(invJ_cache, J_T_cache, detJ_cache)
    end if

    elements=element_count(X)
    
    allocate(invJ_cache(X%dim,X%dim,elements), &
         J_T_cache(X%dim,X%dim,elements), &
         detJ_cache(elements))
#ifdef HAVE_MEMORY_STATS
    call register_allocation("transform_cache", &
         "real", size(invJ_cache)+size(J_T_cache)+size(detJ_cache))
#endif

    x_shape=>ele_shape(X,1)
       
    do ele=1,elements
       X_val=ele_val(X, ele)
       !     |- dx  dx  dx  -|
       !     |  dL1 dL2 dL3  |
       !     |               |
       !     |  dy  dy  dy   |
       ! J = |  dL1 dL2 dL3  |
       !     |               |
       !     |  dz  dz  dz   |
       !     |- dL1 dL2 dL3 -|
       
       ! Form Jacobian.
       ! Since X is linear we need only do this at quadrature point 1.
       J_T_cache(:,:,ele)=matmul(X_val(:,:), x_shape%dn(:, 1, :))
          
       select case (X%dim)
       case(1)
          invJ_cache(:,:,ele)=1.0
       case(2)
          invJ_cache(:,:,ele)=reshape(&
               (/ J_T_cache(2,2,ele),-J_T_cache(1,2,ele),&
               & -J_T_cache(2,1,ele), J_T_cache(1,1,ele)/),(/2,2/))
       case(3)
          ! Calculate (scaled) inverse using recursive determinants.
          forall (i=1:3,k=1:3) 
             invJ_cache(i, k, ele)= &
                  J_T_cache(cyc3(i+1),cyc3(k+1),ele)&
                  &          *J_T_cache(cyc3(i+2),cyc3(k+2),ele) &
                  -J_T_cache(cyc3(i+2),cyc3(k+1),ele)&
                  &          *J_T_cache(cyc3(i+1),cyc3(k+2),ele)
          end forall
       case default
          FLAbort("Unsupported dimension specified.  Universe is 3 dimensional (sorry Albert).")   
       end select
       
       ! Form determinant by expanding minors.
       detJ_cache(ele)=dot_product(J_T_cache(:,1,ele),invJ_cache(:,1,ele))
       
       ! Scale inverse by determinant.
       invJ_cache(:,:,ele)=invJ_cache(:,:,ele)/detJ_cache(ele)
   
    end do
    
  end subroutine construct_cache

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
  
  subroutine transform_to_physical_detwei(X, ele, detwei)
    !!< Fast version of transform_to_physical that only calculates detwei
      
    !! Coordinate field
    type(vector_field), intent(in) :: X
    !! Current element
    integer :: ele
    !! Quadrature weights for physical coordinates.
    real, dimension(:), intent(out):: detwei(:)

    !! Shape function used for coordinate interpolation
    type(element_type), pointer :: x_shape
    !! Column n of X is the position of the nth node. (dim x x_shape%loc)
    !! only need position of n nodes since Jacobian is only calculated once
    real, dimension(X%dim,ele_loc(X,ele)) :: X_val

    real :: J(X%dim, mesh_dim(X)), det
    integer :: gi, dim, ldim
    logical :: x_nonlinear, cache_valid

    x_shape=>ele_shape(X, ele)

    ! Optimisation checks. Optimisations apply to linear elements.
    x_nonlinear= .not.(x_shape%degree==1 .and. x_shape%numbering%family==FAMILY_SIMPLEX)
    
    dim=X%dim ! dimension of space (n/o real coordinates)
    ldim=size(x_shape%dn,3) ! dimension of element (n/o local coordinates)
    if (dim==ldim) then
       
       if ((.not.x_nonlinear).and.cache_transform_elements) then
          cache_valid=retrieve_cached_transform(X, ele, det)
          
          if (cache_valid) then
             detwei=abs(det)*x_shape%quadrature%weight
             return
          end if
          
       end if

!#ifdef DDEBUG
!       if (ele==1) then
!          ewrite(2,*) "Element geometry cache not used."
!       end if
!#endif

       X_val=ele_val(X, ele)
       
       select case (dim)
       case (1)
         do gi=1, x_shape%ngi
           J(1,1)=dot_product(X_val(1,:), x_shape%dn(:,gi,1))
           detwei(gi)=abs(J(1,1))*x_shape%quadrature%weight(gi)
         end do
       case (2)
         do gi=1, x_shape%ngi
            if (x_nonlinear.or.gi==1) then
               ! the Jacobian is the transpose of this
               J=matmul(X_val(:,:), x_shape%dn(:, gi, :))
               ! but that doesn't matter for determinant:
               det=abs(J(1,1)*J(2,2)-J(1,2)*J(2,1))
            end if
           detwei(gi)=det*x_shape%quadrature%weight(gi)
         end do
       case (3)
         do gi=1, x_shape%ngi
            if (x_nonlinear.or.gi==1) then
               ! the Jacobian is the transpose of this
               J=matmul(X_val(:,:), x_shape%dn(:, gi, :))
               ! but that doesn't matter for determinant:
               det=abs( &
                    J(1,1)*(J(2,2)*J(3,3)-J(2,3)*J(3,2)) &
                    -J(1,2)*(J(2,1)*J(3,3)-J(2,3)*J(3,1)) &
                    +J(1,3)*(J(2,1)*J(3,2)-J(2,2)*J(3,1)) &
                    )
            end if
            detwei(gi)= det *x_shape%quadrature%weight(gi)
         end do
       case default
          FLAbort("Unsupported dimension specified.  Universe is 3 dimensional (sorry Albert).")   
       end select
    else if (ldim<dim) then

       X_val=ele_val(X, ele)
      
       ! lower dimensional element (ldim) embedded in higher dimensional space (dim)
       select case (ldim)
       case (1)
          ! 1-dim element embedded in 'dim'-dimensional space:
          do gi=1, x_shape%ngi
             if (x_nonlinear.or.gi==1) then
                ! J is 'dim'-dimensional vector:
                J(:,1)=matmul(X_val(:,:), x_shape%dn(:,gi,1))
                ! length of that
                det=norm2(J(:,1))
             end if
             ! length of that times quad. weight
             detwei(gi)=det*x_shape%quadrature%weight(gi)
          end do
       case (2)
          ! 2-dim element embedded in 'dim'-dimensional space:
          do gi=1, x_shape%ngi
             if (x_nonlinear.or.gi==1) then
                ! J is 2 columns of 2 'dim'-dimensional vectors:
                J=matmul(X_val(:,:), x_shape%dn(:,gi,:))
                ! outer product
                det=abs( &
                     J(2,1)*J(3,2)-J(3,1)*J(2,2) &
                     -J(3,1)*J(1,2)+J(1,1)*J(3,2) &
                     +J(1,1)*J(2,2)-J(2,1)*J(1,2))
             end if
             ! outer product times quad. weight
             detwei(gi)=det *x_shape%quadrature%weight(gi)
          end do
       end select
       
    else
       FLAbort("Don't know how to compute higher-dimensional elements in a lower-dimensional space.")
       
    end if
    
  end subroutine transform_to_physical_detwei
  
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
