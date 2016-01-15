#include "fdebug.h"

subroutine benchmark_parallel_p1_inner_product

  use iso_c_binding, only : c_int8_t

  use libsupermesh_fields, only : triangle_area
  use libsupermesh_fldebug
  use libsupermesh_integer_hash_table
  use libsupermesh_integer_set
  use libsupermesh_halo_ownership
  use libsupermesh_parallel_supermesh
  use libsupermesh_read_halos
  use libsupermesh_read_triangle

  implicit none

#include <finclude/petsc.h90>

  character(len = 1) :: rank_chr
  integer :: ierr, nprocs, rank
  integer :: integer_extent, real_extent
  
  character(len = *), parameter :: basename_a = "data/triangle_0_05_4"
  character(len = *), parameter :: basename_b = "data/square_0_05_4"
 
  integer :: loc_a, nelements_a, nelements_b, nnodes_a, nnodes_b
  integer, parameter :: dim_a = 2, dim_b = 2
  integer, dimension(:), allocatable :: ele_owner_a, ele_owner_b
  integer, dimension(:, :), allocatable, target :: enlist_a
  integer, dimension(:, :), allocatable :: enlist_b
  real, dimension(:), allocatable, target :: field_a
  real, dimension(:), allocatable :: field_b
  real, dimension(:, :), allocatable, target :: positions_a
  real, dimension(:, :), allocatable :: positions_b
  type(halo_type) :: halo_a, halo_b
  
  integer :: data_nelements_a, data_nnodes_a
  integer, dimension(:, :), allocatable, target :: data_enlist_a
  real, dimension(:), allocatable, target :: data_field_a
  real, dimension(:, :), allocatable, target :: data_positions_a
  
  ! Quadrature rule from D. A. Dunavant, "High degree efficient symmetrical
  ! Gaussian quadrature rules for the triangle", International Journal for
  ! Numerical Methods in Engineering, 21, pp. 1129--1148, 1985, appendix II
  real, dimension(3), parameter :: quad_weights = (/1.0D0, 1.0D0, 1.0D0/) / 6.0D0
  real, dimension(3, 3), parameter :: quad_points = reshape((/4.0D0, 1.0D0, 1.0D0, 1.0D0, 4.0D0, 1.0D0, 1.0D0, 1.0D0, 4.0D0/) / 6.0D0, (/3, 3/))
  real :: area_parallel, integral_parallel, real_buffer

  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr);  CHKERRQ(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr);  CHKERRQ(ierr)
  if(nprocs /= 4) then
    FLAbort("Require 4 processes")
  end if    
  write(rank_chr, "(i0)") rank
  call MPI_Type_extent(MPI_INTEGER, integer_extent, ierr);  CHKERRQ(ierr)
  call MPI_Type_extent(MPI_DOUBLE_PRECISION, real_extent, ierr);  CHKERRQ(ierr)
  
  call read_node(basename_a // "_" // rank_chr // ".node", dim = dim_a, coords = positions_a)
  call read_ele(basename_a // "_" // rank_chr // ".ele", dim = dim_a, enlist = enlist_a)
  nnodes_a = size(positions_a, 2)
  loc_a = size(enlist_a, 1)
  nelements_a = size(enlist_a, 2)
  call read_halo(basename_a, halo_a, level = 2)
  allocate(ele_owner_a(nelements_a))
  call element_ownership(nnodes_a, enlist_a, halo_a, ele_owner_a)
  allocate(field_a(nnodes_a))
  field_a = positions_a(1, :)
  
  call read_node(basename_b // "_"// rank_chr // ".node", dim = dim_b, coords = positions_b)
  call read_ele(basename_b // "_" // rank_chr // ".ele", dim = dim_b, enlist = enlist_b)
  nnodes_b = size(positions_b, 2)
  nelements_b = size(enlist_b, 2)
  call read_halo(basename_b, halo_b, level = 2)
  allocate(ele_owner_b(nelements_b))
  call element_ownership(nnodes_b, enlist_b, halo_b, ele_owner_b)
  allocate(field_b(nnodes_b))
  field_b = positions_b(2, :)
  
  area_parallel = 0.0D0
  integral_parallel = 0.0D0
  call parallel_supermesh(positions_b, enlist_b, ele_owner_b, &
                        & positions_a, enlist_a, ele_owner_a, &
                        & donor_ele_data, unpack_data_a, intersection_calculation, &
                        & comm = MPI_COMM_WORLD)
  call cleanup_unpack_data_a()
  call MPI_Allreduce(area_parallel, real_buffer, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  area_parallel = real_buffer
  call MPI_Allreduce(integral_parallel, real_buffer, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  integral_parallel = real_buffer
  if(rank == 0) then
    print *, "Area     = ", area_parallel
    print *, "Integral = ", integral_parallel
  end if
                        
  deallocate(positions_a, enlist_a, positions_b, enlist_b, ele_owner_a, ele_owner_b, field_a, field_b)
  call deallocate(halo_a)
  call deallocate(halo_b)

contains

  subroutine donor_ele_data(eles_a, data_a)
    integer, dimension(:), intent(in)                   :: eles_a
    integer(kind = c_int8_t), dimension(:), allocatable :: data_a
    
    integer :: data_nnodes_a, i, lnode, ndata_a, node, position
    integer, dimension(:, :), allocatable :: data_enlist_a
    real, dimension(:), allocatable :: data_field_a
    real, dimension(:, :), allocatable :: data_positions_a
    type(integer_hash_table) :: node_map
    type(integer_set) :: nodes
    
    call allocate(nodes)
    do i = 1, size(eles_a)
      do lnode = 1, loc_a
        call insert(nodes, enlist_a(lnode, eles_a(i)))
      end do
    end do
    data_nnodes_a = key_count(nodes)
    call allocate(node_map)
    allocate(data_positions_a(dim_a, data_nnodes_a), data_field_a(data_nnodes_a))
    do i = 1, data_nnodes_a
      node = fetch(nodes, i)
      call insert(node_map, node, i)
      data_positions_a(:, i) = positions_a(:, node)
      data_field_a(i) = field_a(node)
    end do
    call deallocate(nodes)
    allocate(data_enlist_a(loc_a, size(eles_a)))
    do i = 1, size(eles_a)
      do lnode = 1, loc_a
        data_enlist_a(lnode, i) = fetch(node_map, enlist_a(lnode, eles_a(i)))
      end do
    end do
    call deallocate(node_map)
    
    ndata_a = (2 + loc_a * size(eles_a)) * integer_extent + (dim_a + 1) * size(data_field_a) * real_extent
    allocate(data_a(ndata_a))
    position = 0
    call MPI_Pack((/size(eles_a), size(data_field_a)/), 2, MPI_INTEGER, data_a, ndata_a, position, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    call MPI_Pack(data_enlist_a, loc_a * size(eles_a), MPI_INTEGER, data_a, ndata_a, position, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    call MPI_Pack(data_positions_a, dim_a * size(data_field_a), MPI_DOUBLE_PRECISION, data_a, ndata_a, position, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    call MPI_Pack(data_field_a, size(data_field_a), MPI_DOUBLE_PRECISION, data_a, ndata_a, position, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    
    deallocate(data_enlist_a, data_positions_a, data_field_a)
    
  end subroutine donor_ele_data
  
  subroutine unpack_data_a(data_a)
    integer(kind = c_int8_t), dimension(:), intent(in) :: data_a
    
    integer :: position
    integer, dimension(2) :: data_header_a
    
    call cleanup_unpack_data_a()
    
    position = 0
    call MPI_Unpack(data_a, size(data_a), position, data_header_a, 2, MPI_INTEGER, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    data_nelements_a = data_header_a(1)
    data_nnodes_a = data_header_a(2)
    allocate(data_enlist_a(loc_a, data_nelements_a), data_positions_a(dim_a, data_nnodes_a), data_field_a(data_nnodes_a))
    call MPI_Unpack(data_a, size(data_a), position, data_enlist_a, loc_a * data_nelements_a, MPI_INTEGER, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    call MPI_Unpack(data_a, size(data_a), position, data_positions_a, dim_a * data_nnodes_a, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    call MPI_Unpack(data_a, size(data_a), position, data_field_a, data_nnodes_a, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    
  end subroutine unpack_data_a
  
  subroutine cleanup_unpack_data_a()
    if(allocated(data_enlist_a)) then
      deallocate(data_enlist_a, data_positions_a, data_field_a)
    end if
  end subroutine cleanup_unpack_data_a
  
  pure function basis_functions(cell_coords, coord) result(fns)
    real, dimension(2, 3), intent(in) :: cell_coords
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
    fns(1) = 1.0 - fns(2) - fns(3)

  end function basis_functions

  pure function interpolate(cell_coords_d, cell_x_d, coord_s) result(x_s)
    real, dimension(2, 3), intent(in) :: cell_coords_d
    real, dimension(3), intent(in) :: cell_x_d
    real, dimension(2), intent(in) :: coord_s

    real :: x_s

    x_s = dot_product(basis_functions(cell_coords_d, coord_s), cell_x_d)

  end function interpolate
  
  subroutine intersection_calculation(positions_c, ele_b, ele_a, local)
    ! dim x loc_c x nelements_c
    real, dimension(:, :, :), intent(in) :: positions_c
    integer, intent(in) :: ele_b
    integer, intent(in) :: ele_a
    logical, intent(in) :: local
    
    integer :: ele_c, i, nelements_c
    integer, dimension(:, :), pointer :: lenlist_a
    real :: area
    real, dimension(2) :: quad_point
    real, dimension(:), pointer :: lfield_a
    real, dimension(:, :), pointer :: lpositions_a
    
    nelements_c = size(positions_c, 3)
    
    if(local) then
      lenlist_a => enlist_a
      lpositions_a => positions_a
      lfield_a => field_a
    else
      lenlist_a => data_enlist_a
      lpositions_a => data_positions_a
      lfield_a => data_field_a
    end if
    
    do ele_c = 1, nelements_c
      area = triangle_area(positions_c(:, :, ele_c))
      area_parallel = area_parallel + area
      do i = 1, size(quad_weights)
        quad_point(1) = dot_product(quad_points(:, i), positions_c(1, :, ele_c))
        quad_point(2) = dot_product(quad_points(:, i), positions_c(2, :, ele_c))
        integral_parallel = integral_parallel + 2.0D0 * quad_weights(i) * area &
                                              & * interpolate(lpositions_a(:, lenlist_a(:, ele_a)), lfield_a(lenlist_a(:, ele_a)), quad_point) &
                                              & * interpolate(positions_b(:, enlist_b(:, ele_b)), field_b(enlist_b(:, ele_b)), quad_point)
      end do
    end do
    
  end subroutine intersection_calculation
  
end subroutine benchmark_parallel_p1_inner_product
