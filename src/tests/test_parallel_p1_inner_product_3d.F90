#include "libsupermesh_debug.h"

subroutine test_parallel_p1_inner_product_3d() bind(c)

  use iso_c_binding, only : c_int8_t
  use iso_fortran_env, only : error_unit, output_unit

  use libsupermesh_debug, only : abort_pinpoint
  use libsupermesh_halo_ownership, only : element_ownership
  use libsupermesh_parallel_supermesh, only : parallel_supermesh
  use libsupermesh_read_halos, only : halo_type, deallocate, read_halo
  use libsupermesh_read_triangle, only : read_ele, read_node
  use libsupermesh_supermesh, only : tetrahedron_volume
  use libsupermesh_unittest_tools, only : operator(.fne.), report_test

  implicit none

#include <mpif.h>
  
  ! Input Triangle mesh base names
  character(len = *), parameter :: basename_a = "data/pyramid_0_05", &
                                 & basename_b = "data/cube_0_05"
  real, parameter :: volume_ref = 3.3333333333333331D+02 / 1.0D3, integral_ref = 1.2499999999999989D+04 / 1.0D5

  character(len = int(log10(real(huge(0)))) + 1) :: rank_chr, nprocs_chr
  integer :: ierr, integer_extent, nprocs, rank, real_extent
 
  integer :: nelements_a, nelements_b, nnodes_a, nnodes_b
  integer, dimension(:), allocatable :: ele_owner_a, ele_owner_b
  integer, dimension(:, :), allocatable :: enlist_a, enlist_b
  real, dimension(:), allocatable :: field_a
  real, dimension(:), allocatable, target :: field_b
  real, dimension(:, :), allocatable :: positions_a, positions_b
  type(halo_type) :: halo_a, halo_b
  
  integer :: data_nnodes_b
  real, dimension(:), allocatable, target :: data_field_b
  
  ! Quadrature rule from "Approximate Calculation of Multiple Integrals",
  ! A. H. Stroud, Prentice-Hall, Inc., 1971, section 8.8, T_n: 2-1 upper sign.
  ! See also P. C. Hammer and A. H. Stroud, "Numerical integration over
  ! simplexes", Mathematical Tables and Other Aids to Computation, 10,
  ! pp. 137--139, 1956.
  real, dimension(4), parameter :: quad_weights = (/1.0D0, 1.0D0, 1.0D0, 1.0D0/) / 4.0D0
  real, dimension(4, 4), parameter :: quad_points = reshape((/5.8541019662496852D-01, 1.3819660112501050D-01, 1.3819660112501050D-01, 1.3819660112501050D-01, &
                                                            & 1.3819660112501050D-01, 5.8541019662496852D-01, 1.3819660112501050D-01, 1.3819660112501050D-01, &
                                                            & 1.3819660112501050D-01, 1.3819660112501050D-01, 5.8541019662496852D-01, 1.3819660112501050D-01, &
                                                            & 1.3819660112501050D-01, 1.3819660112501050D-01, 1.3819660112501050D-01, 5.8541019662496852D-01/), (/4, 4/))
  real :: volume_parallel, integral_parallel

  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr);  assert(ierr == MPI_SUCCESS)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr);  assert(ierr == MPI_SUCCESS)
  write(rank_chr, "(i0)") rank
  write(nprocs_chr, "(i0)") nprocs
  rank_chr = adjustl(rank_chr)
  nprocs_chr = adjustl(nprocs_chr)
  call MPI_Type_extent(MPI_INTEGER, integer_extent, ierr);  assert(ierr == MPI_SUCCESS)
  call MPI_Type_extent(MPI_DOUBLE_PRECISION, real_extent, ierr);  assert(ierr == MPI_SUCCESS)

  ! Read the donor mesh partition
  call read_node(trim(basename_a) // "_" // trim(nprocs_chr) // "_" // trim(rank_chr) // ".node", dim = 3, coords = positions_a)
  call read_ele(trim(basename_a) // "_" // trim(nprocs_chr) // "_" // trim(rank_chr) // ".ele", dim = 3, enlist = enlist_a)
  nnodes_a = size(positions_a, 2)
  nelements_a = size(enlist_a, 2)
  ! Read donor mesh halo data ...
  call read_halo(trim(basename_a) // "_" // trim(nprocs_chr), halo_a, level = 2)
  ! ... and determine the donor mesh element ownership
  allocate(ele_owner_a(nelements_a))
  call element_ownership(nnodes_a, enlist_a, halo_a, ele_owner_a)
  ! Construct a donor P1 field equal to: x + z
  allocate(field_a(nnodes_a))
  field_a = positions_a(1, :) + positions_a(3, :)

  ! Read the target mesh partition
  call read_node(trim(basename_b) // "_" // trim(nprocs_chr) // "_" // trim(rank_chr) // ".node", dim = 3, coords = positions_b)
  call read_ele(trim(basename_b) // "_" // trim(nprocs_chr) // "_" // trim(rank_chr) // ".ele", dim = 3, enlist = enlist_b)
  nnodes_b = size(positions_b, 2)
  nelements_b = size(enlist_b, 2)
  ! Read target mesh halo data ...
  call read_halo(trim(basename_b) // "_" // trim(nprocs_chr), halo_b, level = 2)
  ! ... and determine the target mesh element ownership
  allocate(ele_owner_b(nelements_b))
  call element_ownership(nnodes_b, enlist_b, halo_b, ele_owner_b)
  ! Construct a target P1 field equal to: y
  allocate(field_b(nnodes_b))
  field_b = positions_b(2, :)

  ! Perform multi-mesh integration
  volume_parallel = 0.0D0
  integral_parallel = 0.0D0
  call parallel_supermesh(positions_a, enlist_a, ele_owner_a, &
                        & positions_b, enlist_b, ele_owner_b, &
                        & pack_data_b, unpack_data_b, intersection_calculation, &
                        & comm = MPI_COMM_WORLD)
  ! Deallocate any remaining unpacked communicated data
  call cleanup_unpack_data_b()

  ! Sum all process contributions to the multi-mesh integrals
  call MPI_Allreduce(MPI_IN_PLACE, volume_parallel, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
  call MPI_Allreduce(MPI_IN_PLACE, integral_parallel, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)

  flush(output_unit)
  flush(error_unit)
  call MPI_Barrier(MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
  if(rank == 0) then
    ! Display the multi-mesh integrals on rank 0
    print "(a,e26.18e3,a,e26.18e3,a)", "Volume   = ", volume_parallel, " (error = ", abs(volume_parallel - volume_ref), ")"
    print "(a,e26.18e3,a,e26.18e3,a)", "Integral = ", integral_parallel, " (error = ", abs(integral_parallel - integral_ref), ")"
    
    ! Test the multi-mesh integrals against the reference values
    call report_test("[test_parallel_p1_inner_product_3d volume]", volume_parallel .fne. volume_ref, .false., "Incorrect area")
    call report_test("[test_parallel_p1_inner_product_3d integral]", integral_parallel .fne. integral_ref, .false., "Incorrect integral")
  end if
  flush(output_unit)
  flush(error_unit)

  ! Cleanup
  deallocate(positions_a, enlist_a, ele_owner_a, field_a, &
           & positions_b, enlist_b, ele_owner_b, field_b)
  call deallocate(halo_a)
  call deallocate(halo_b)

contains

  ! Given the provided mesh vertices and elements, pack data for communication
  subroutine pack_data_b(nodes_b, eles_b, data_b)
    ! Mesh vertices to be communicated
    integer, dimension(:), intent(in) :: nodes_b
    ! Mesh elements to be communicated
    integer, dimension(:), intent(in) :: eles_b
    ! Packed data for communication
    integer(kind = c_int8_t), dimension(:), allocatable, intent(out) :: data_b
    
    integer :: ndata_b, position
    real, dimension(:), allocatable :: data_field_b
    
    ! Gather P1 field values for communication
    allocate(data_field_b(size(nodes_b)))
    data_field_b = field_b(nodes_b)
    
    ! Pack data for communication:
    !   (number of P1 nodes) reals -- communicated P1 field values
    ndata_b = size(data_field_b) * real_extent
    allocate(data_b(ndata_b))
    position = 0
    call MPI_Pack(data_field_b, size(data_field_b), MPI_DOUBLE_PRECISION, data_b, ndata_b, position, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
    
    deallocate(data_field_b)
    
  end subroutine pack_data_b
  
  ! Unpack communicated data
  subroutine unpack_data_b(nnodes_b, nelements_b, data_b)
    ! Number of communicated mesh vertices
    integer, intent(in) :: nnodes_b
    ! Number of communicated elements
    integer, intent(in) :: nelements_b
    integer(kind = c_int8_t), dimension(:), intent(in) :: data_b
    
    integer :: position
    
    ! Deallocate any previously unpacked communicated data
    call cleanup_unpack_data_b()
    
    position = 0
    ! Store the number of P1 nodes
    data_nnodes_b = nnodes_b
    ! Unpack the P1 field values
    allocate(data_field_b(data_nnodes_b))
    call MPI_Unpack(data_b, size(data_b), position, data_field_b, data_nnodes_b, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr);  assert(ierr == MPI_SUCCESS)
    
  end subroutine unpack_data_b
  
  ! Deallocate any previously unpacked communicated data
  subroutine cleanup_unpack_data_b()
    if(allocated(data_field_b)) deallocate(data_field_b)
  end subroutine cleanup_unpack_data_b
  
  ! Evaluate P1 basis functions at a given point
  pure function basis_functions_p1(cell_coords, coord) result(fns)
    ! Tetrahedron vertex coordinates
    ! Shape: dim x loc_p1
    real, dimension(3, 4), intent(in) :: cell_coords
    ! Coordinate at which to evaluate the basis functions
    ! Shape: dim
    real, dimension(3), intent(in) :: coord

    real, dimension(4) :: fns

    integer :: i, j
    real, dimension(4) :: tmp
    real, dimension(3, 4) :: mat
    
    mat(:, 1:3) = cell_coords(:, 2:4) - spread(cell_coords(:, 1), 2, 3)
    mat(:, 4) = coord - cell_coords(:, 1)
    ! GE with partial pivoting
    do i = 1, size(mat, 1)
      j = maxloc(abs(mat(i:, i)), 1) + i - 1
      if(j /= i) then
        tmp = mat(i, :)
        mat(i, :) = mat(j, :)
        mat(j, :) = tmp
      end if
      mat(i, :) = mat(i, :) / mat(i, i)
      do j = i + 1, size(mat, 1)
        mat(j, :) = mat(j, :) - mat(j, i) * mat(i, :)
      end do
    end do
    do i = size(mat, 1), 1, -1
      do j = i - 1, 1, -1
        mat(j, :) = mat(j, :) - mat(j, i) * mat(i, :)
      end do
    end do
    
    fns(2:4) = mat(:, 4)
    fns(1) = 1.0D0 - fns(2) - fns(3) - fns(4)

  end function basis_functions_p1

  ! Interpolate a P1 function at given point
  pure function interpolate_p1(cell_coords_d, cell_x_d, coord_s) result(x_s)
    ! Tetrahedron vertex coordinates
    ! Shape: dim x loc_p1
    real, dimension(3, 4), intent(in) :: cell_coords_d
    ! P1 nodal values
    ! Shape: loc_p1
    real, dimension(4), intent(in) :: cell_x_d
    ! Coordinate at which to evaluate the P1 function
    ! Shape: dim
    real, dimension(3), intent(in) :: coord_s

    real :: x_s

    x_s = dot_product(basis_functions_p1(cell_coords_d, coord_s), cell_x_d)

  end function interpolate_p1
  
  ! Perform calculations on the local supermesh
  subroutine intersection_calculation(positions_a, positions_b, positions_c, nodes_b, ele_a, ele_b, local)
    ! Target mesh element vertex coordinates
    ! Shape: dim x loc_a
    real, dimension(:, :), intent(in) :: positions_a
    ! Donor mesh element vertex coordinates
    ! Shape: dim x loc_b
    real, dimension(:, :), intent(in) :: positions_b
    ! Supermesh element vertex coordinates
    ! Shape: dim x loc_c x nelements_c
    real, dimension(:, :, :), intent(in) :: positions_c
    ! Donor mesh vertex indices
    ! Shape: loc_b
    integer, dimension(:), intent(in) :: nodes_b
    ! Target mesh element
    integer, intent(in) :: ele_a
    ! Donor mesh element
    integer, intent(in) :: ele_b
    ! Whether this is a local calculation or a calculation using communicated
    ! data
    logical, intent(in) :: local
    
    integer :: ele_c, i
    real :: volume
    real, dimension(3) :: quad_point
    real, dimension(:), pointer :: lfield_b
    
    if(local) then
      ! If this is a local calculation, use the local P1 field data
      lfield_b => field_b
    else
      ! Otherwise, use the unpacked communicated P1 field data
      lfield_b => data_field_b
    end if
    
    do ele_c = 1, size(positions_c, 3)
      ! Compute the supermesh tetrahedron volume
      volume = tetrahedron_volume(positions_c(:, :, ele_c))
      ! Local contribution to the intersection volume
      volume_parallel = volume_parallel + volume
      ! Local contribution to the multi-mesh inner product, evaluated using
      ! degree 2 quadrature
      do i = 1, size(quad_weights)
        quad_point(1) = dot_product(quad_points(:, i), positions_c(1, :, ele_c))
        quad_point(2) = dot_product(quad_points(:, i), positions_c(2, :, ele_c))
        quad_point(3) = dot_product(quad_points(:, i), positions_c(3, :, ele_c))
        integral_parallel = integral_parallel + quad_weights(i) * volume &
                                              & * interpolate_p1(positions_b, lfield_b(nodes_b), quad_point) &
                                              & * interpolate_p1(positions_a, field_a(enlist_a(:, ele_a)), quad_point)
      end do
    end do
        
  end subroutine intersection_calculation
  
end subroutine test_parallel_p1_inner_product_3d