subroutine benchmark_parallel_p1_inner_product

  use iso_c_binding, only : c_int8_t

  use libsupermesh_fields, only : triangle_area
  use libsupermesh_halo_ownership, only : element_ownership
  use libsupermesh_parallel_supermesh, only : parallel_supermesh
  use libsupermesh_read_halos, only : halo_type, deallocate, read_halo
  use libsupermesh_read_triangle, only : read_ele, read_node

  implicit none

#include <finclude/petsc.h90>
  
  ! Input Triangle mesh base names
  character(len = *), parameter :: basename_a = "data/triangle_0_01_4", &
                                 & basename_b = "data/square_0_01_4"

  character(len = int(log10(real(huge(0)))) + 2) :: rank_chr
  integer :: ierr, integer_extent, rank, real_extent
  real :: real_buffer
 
  integer :: nelements_a, nelements_b, nnodes_a, nnodes_b
  integer, dimension(:), allocatable :: ele_owner_a, ele_owner_b
  integer, dimension(:, :), allocatable :: enlist_a, enlist_b
  real, dimension(:), allocatable :: field_b
  real, dimension(:), allocatable, target :: field_a
  real, dimension(:, :), allocatable :: positions_a, positions_b
  type(halo_type) :: halo_a, halo_b
  
  integer :: data_nnodes_a
  real, dimension(:), allocatable, target :: data_field_a
  
  ! Quadrature rule from D. A. Dunavant, "High degree efficient symmetrical
  ! Gaussian quadrature rules for the triangle", International Journal for
  ! Numerical Methods in Engineering, 21, pp. 1129--1148, 1985, appendix II
  real, dimension(3), parameter :: quad_weights = (/1.0D0, 1.0D0, 1.0D0/) / 3.0D0
  real, dimension(3, 3), parameter :: quad_points = reshape((/4.0D0, 1.0D0, 1.0D0, 1.0D0, 4.0D0, 1.0D0, 1.0D0, 1.0D0, 4.0D0/) / 6.0D0, (/3, 3/))
  real :: area_parallel, integral_parallel

  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr);  CHKERRQ(ierr)
  write(rank_chr, "(i0)") rank
  rank_chr = adjustl(rank_chr)
  call MPI_Type_extent(MPI_INTEGER, integer_extent, ierr);  CHKERRQ(ierr)
  call MPI_Type_extent(MPI_DOUBLE_PRECISION, real_extent, ierr);  CHKERRQ(ierr)
  
  ! Read the donor mesh partition
  call read_node(trim(basename_a) // "_" // trim(rank_chr) // ".node", dim = 2, coords = positions_a)
  call read_ele(trim(basename_a) // "_" // trim(rank_chr) // ".ele", dim = 2, enlist = enlist_a)
  nnodes_a = size(positions_a, 2)
  nelements_a = size(enlist_a, 2)
  ! Read donor mesh halo data ...
  call read_halo(basename_a, halo_a, level = 2)
  ! ... and determine the donor mesh element ownership
  allocate(ele_owner_a(nelements_a))
  call element_ownership(nnodes_a, enlist_a, halo_a, ele_owner_a)
  ! Construct a donor P1 field equal to: x
  allocate(field_a(nnodes_a))
  field_a = positions_a(1, :)
  
  ! Read the target mesh partition
  call read_node(trim(basename_b) // "_"// trim(rank_chr) // ".node", dim = 2, coords = positions_b)
  call read_ele(trim(basename_b) // "_" // trim(rank_chr) // ".ele", dim = 2, enlist = enlist_b)
  nnodes_b = size(positions_b, 2)
  nelements_b = size(enlist_b, 2)
  ! Read target mesh halo data ...
  call read_halo(basename_b, halo_b, level = 2)
  ! ... and determine the target mesh element ownership
  allocate(ele_owner_b(nelements_b))
  call element_ownership(nnodes_b, enlist_b, halo_b, ele_owner_b)
  ! Construct a target P1 field equal to: y
  allocate(field_b(nnodes_b))
  field_b = positions_b(2, :)
  
  ! Perform multi-mesh integration
  area_parallel = 0.0D0
  integral_parallel = 0.0D0
  call parallel_supermesh(positions_b, enlist_b, ele_owner_b, &
                        & positions_a, enlist_a, ele_owner_a, &
                        & donor_ele_data, unpack_data_a, intersection_calculation, &
                        & comm = MPI_COMM_WORLD)
  ! Deallocate any remaining unpacked communicated data
  call cleanup_unpack_data_a()
  ! Sum all process contributions to the multi-mesh integrals
  call MPI_Allreduce(area_parallel, real_buffer, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  area_parallel = real_buffer
  call MPI_Allreduce(integral_parallel, real_buffer, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  integral_parallel = real_buffer
  if(rank == 0) then
    ! Display the multi-mesh integrals on rank 0
    print "(a,e26.18e3,a,e26.18e3,a)", "Area     = ", area_parallel, " (error = ", abs(area_parallel - 0.5D0), ")"
    print "(a,e26.18e3,a,e26.18e3,a)", "Integral = ", integral_parallel, " (error = ", abs(integral_parallel - 8.3333333333333398D-02), ")"
  end if
                        
  deallocate(positions_a, enlist_a, ele_owner_a, field_a, &
           & positions_b, enlist_b, ele_owner_b, field_b)
  call deallocate(halo_a)
  call deallocate(halo_b)

contains

  ! Given the provided mesh vertices and elements, pack data for communication
  subroutine donor_ele_data(nodes_a, eles_a, data_a)
    ! Mesh vertices to be communicated
    integer, dimension(:), intent(in) :: nodes_a
    ! Mesh elements to be communicated
    integer, dimension(:), intent(in) :: eles_a
    ! Packed data for communication
    integer(kind = c_int8_t), dimension(:), allocatable, intent(out) :: data_a
    
    integer :: ndata_a, position
    real, dimension(:), allocatable :: data_field_a
    
    ! Gather P1 field values for communication
    allocate(data_field_a(size(nodes_a)))
    data_field_a = field_a(nodes_a)
    
    ! Pack data for communication:
    !   1 integer                  -- number of P1 nodes
    !   (number of P1 nodes) reals -- communicated P1 field values
    ndata_a = integer_extent + size(data_field_a) * real_extent
    allocate(data_a(ndata_a))
    position = 0
    call MPI_Pack(size(data_field_a), 1, MPI_INTEGER, data_a, ndata_a, position, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    call MPI_Pack(data_field_a, size(data_field_a), MPI_DOUBLE_PRECISION, data_a, ndata_a, position, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    
    deallocate(data_field_a)
    
  end subroutine donor_ele_data
  
  ! Unpack communicated data
  subroutine unpack_data_a(data_a)
    integer(kind = c_int8_t), dimension(:), intent(in) :: data_a
    
    integer :: position
    
    ! Deallocate any previously unpacked communicated data
    call cleanup_unpack_data_a()
    
    position = 0
    ! Unpack the number of P1 nodes
    call MPI_Unpack(data_a, size(data_a), position, data_nnodes_a, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    ! Unpack the P1 field values
    allocate(data_field_a(data_nnodes_a))
    call MPI_Unpack(data_a, size(data_a), position, data_field_a, data_nnodes_a, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    
  end subroutine unpack_data_a
  
  ! Deallocate any previously unpacked communicated data
  subroutine cleanup_unpack_data_a()
    if(allocated(data_field_a)) deallocate(data_field_a)
  end subroutine cleanup_unpack_data_a
  
  ! Evaluate P1 basis functions at a given point
  pure function basis_functions_p1(cell_coords, coord) result(fns)
    ! Triangle vertex coordinates
    ! Shape: dim x loc_p1
    real, dimension(2, 3), intent(in) :: cell_coords
    ! Coordinate at which to evaluate the basis functions
    ! Shape: dim
    real, dimension(2), intent(in) :: coord

    real, dimension(3) :: fns

    real, dimension(2) :: e1, e2
    real, dimension(2, 2) :: jac
        
    e1 = cell_coords(:, 2) - cell_coords(:, 1)
    e2 = cell_coords(:, 3) - cell_coords(:, 1)

    jac(1, 1) =  e2(2);  jac(1, 2) = -e2(1)
    jac(2, 1) = -e1(2);  jac(2, 2) =  e1(1)
    jac = jac / (jac(1, 1) * jac(2, 2) - jac(1, 2) * jac(2, 1))

    fns(2:3) = matmul(jac, coord - cell_coords(:, 1))
    fns(1) = 1.0D0 - fns(2) - fns(3)

  end function basis_functions_p1

  ! Interpolate a P1 function at given point
  pure function interpolate_p1(cell_coords_d, cell_x_d, coord_s) result(x_s)
    ! Triangle vertex coordinates
    ! Shape: dim x loc_p1
    real, dimension(2, 3), intent(in) :: cell_coords_d
    ! P1 nodal values
    ! Shape: loc_p1
    real, dimension(3), intent(in) :: cell_x_d
    ! Coordinate at which to evaluate the P1 function
    ! Shape: dim
    real, dimension(2), intent(in) :: coord_s

    real :: x_s

    x_s = dot_product(basis_functions_p1(cell_coords_d, coord_s), cell_x_d)

  end function interpolate_p1
  
  ! Perform calculations on the local supermesh
  subroutine intersection_calculation(positions_b, positions_a, positions_c, nodes_a, ele_b, ele_a, local)
    ! Target mesh element vertex coordinates
    ! Shape: dim x loc_b
    real, dimension(:, :), intent(in) :: positions_b
    ! Donor mesh element vertex coordniates
    ! Shape: dim x loc_a
    real, dimension(:, :), intent(in) :: positions_a
    ! Supermesh element vertex coordinates
    ! Shape: dim x loc_c x nelements_c
    real, dimension(:, :, :), intent(in) :: positions_c
    ! Donor mesh vertex indices
    ! Shape: loc_a
    integer, dimension(:), intent(in) :: nodes_a
    ! Target mesh element
    integer, intent(in) :: ele_b
    ! Donor mesh element
    integer, intent(in) :: ele_a
    ! Whether this is a local calculation or a calculation using communicated
    ! data
    logical, intent(in) :: local
    
    integer :: ele_c, i
    real :: area
    real, dimension(2) :: quad_point
    real, dimension(:), pointer :: lfield_a
    
    if(local) then
      ! If this is a local calculation, use the local P1 field data
      lfield_a => field_a
    else
      ! Otherwise, use the unpacked communicated P1 field data
      lfield_a => data_field_a
    end if
    
    do ele_c = 1, size(positions_c, 3)
      ! Compute the supermesh triangle area
      area = triangle_area(positions_c(:, :, ele_c))
      ! Local contribution to the intersection area
      area_parallel = area_parallel + area
      ! Local contribution to the multi-mesh inner product, evaluated using
      ! degree 2 quadrature
      do i = 1, size(quad_weights)
        quad_point(1) = dot_product(quad_points(:, i), positions_c(1, :, ele_c))
        quad_point(2) = dot_product(quad_points(:, i), positions_c(2, :, ele_c))
        integral_parallel = integral_parallel + quad_weights(i) * area &
                                              & * interpolate_p1(positions_a, lfield_a(nodes_a), quad_point) &
                                              & * interpolate_p1(positions_b, field_b(enlist_b(:, ele_b)), quad_point)
      end do
    end do
        
  end subroutine intersection_calculation
  
end subroutine benchmark_parallel_p1_inner_product
