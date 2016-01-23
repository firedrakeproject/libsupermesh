subroutine benchmark_parallel_p2_inner_product

  use iso_c_binding, only : c_int8_t

  use libsupermesh_fields, only : triangle_area
  use libsupermesh_halo_ownership, only : element_ownership
  use libsupermesh_integer_hash_table, only : integer_hash_table, allocate, &
    & deallocate, has_key, fetch, insert
  use libsupermesh_integer_set, only : integer_set, allocate, deallocate, &
    & insert, key_count, fetch
  use libsupermesh_parallel_supermesh, only : parallel_supermesh, get_times
  use libsupermesh_read_halos, only : halo_type, deallocate, read_halo
  use libsupermesh_read_triangle, only : read_ele, read_node

  implicit none

#include <finclude/petsc.h90>
#include "fdebug.h"
  
  ! Input Triangle mesh base names
  character(len = *), parameter :: basename_a = "data/triangle_0_01", &
                                 & basename_b = "data/square_0_01"

  character(len = int(log10(real(huge(0)))) + 2) :: rank_chr, size_chr
  integer :: ierr, integer_extent, rank, real_extent
  real :: real_buffer

  integer :: nelements_a, nelements_b, nnodes_p1_a, nnodes_p2_a, &
    & nnodes_p1_b, nnodes_p2_b
  integer, dimension(:), allocatable :: ele_owner_a, ele_owner_b
  integer, dimension(:, :), allocatable :: enlist_p1_a, enlist_p1_b, &
    & enlist_p2_b
  integer, dimension(:, :), allocatable, target :: enlist_p2_a
  real, dimension(:), allocatable :: interpolated, field_b
  real, dimension(:), allocatable, target :: field_a
  real, dimension(:, :), allocatable :: positions_a, positions_b
  type(halo_type) :: halo_a, halo_b

  integer :: data_nelements_a, data_nnodes_p2_a
  integer, dimension(:, :), allocatable, target :: data_enlist_p2_a
  real, dimension(:), allocatable, target :: data_field_a

  real :: t0, parallel_read_time, parallel_interpolation_time, parallel_time
  real :: parallel_time_read_min, parallel_time_read_max, parallel_time_read_sum
  real :: parallel_time_min, parallel_time_max, parallel_time_sum
  real :: parallel_time_interpolation_min, parallel_time_interpolation_max, parallel_time_interpolation_sum
  integer :: nprocs
#if PROFILE == 1
  real :: all_to_all_max, all_to_all_min, all_to_all_sum
  real :: point_to_point_max, point_to_point_min, point_to_point_sum
#endif

  ! P2 mass matrix
  real, dimension(6, 6), parameter :: mass_p2 = reshape((/ 6.0D0, -1.0D0, -1.0D0,  0.0D0, -4.0D0,  0.0D0, &
                                                        & -1.0D0,  6.0D0, -1.0D0,  0.0D0,  0.0D0, -4.0D0, &
                                                        & -1.0D0, -1.0D0,  6.0D0, -4.0D0,  0.0D0,  0.0D0, &
                                                        &  0.0D0,  0.0D0, -4.0D0, 32.0D0, 16.0D0, 16.0D0, &
                                                        & -4.0D0,  0.0D0,  0.0D0, 16.0D0, 32.0D0, 16.0D0, &
                                                        &  0.0D0, -4.0D0,  0.0D0, 16.0D0, 16.0D0, 32.0D0/) / 360.0D0, (/6, 6/))
  real :: area_parallel, integral_parallel

  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr);  CHKERRQ(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr); CHKERRQ(ierr)
  write(rank_chr, "(i0)") rank
  write(size_chr, "(i0)") nprocs
  rank_chr = adjustl(rank_chr)
  size_chr = adjustl(size_chr)
  call MPI_Type_extent(MPI_INTEGER, integer_extent, ierr);  CHKERRQ(ierr)
  call MPI_Type_extent(MPI_DOUBLE_PRECISION, real_extent, ierr);  CHKERRQ(ierr)

  t0 = mpi_wtime()
  ! Read the donor mesh partition
  call read_node(trim(basename_a) // "_" // trim(size_chr) // "_" // trim(rank_chr) // ".node", dim = 2, coords = positions_a)
  call read_ele(trim(basename_a) // "_" // trim(size_chr) // "_" // trim(rank_chr) // ".ele", dim = 2, enlist = enlist_p1_a)
  nnodes_p1_a = size(positions_a, 2)
  nelements_a = size(enlist_p1_a, 2)
  ! Read donor mesh halo data ...
  call read_halo(trim(basename_a) // "_" // trim(size_chr), halo_a, level = 2)
  parallel_read_time = mpi_wtime() - t0

  t0 = mpi_wtime()
  ! ... and determine the donor mesh element ownership
  allocate(ele_owner_a(nelements_a))
  call element_ownership(nnodes_p1_a, enlist_p1_a, halo_a, ele_owner_a)
  ! Generate the donor P2 element-node graph
  allocate(enlist_p2_a(6, nelements_a))
  call p2_connectivity(nnodes_p1_a, enlist_p1_a, nnodes_p2_a, enlist_p2_a)
  ! Construct a donor P2 field equal to: x y
  allocate(field_a(nnodes_p2_a), interpolated(nnodes_p2_a))
  call interpolate_p1_p2(enlist_p1_a, positions_a(1, :), enlist_p2_a, field_a)
  call interpolate_p1_p2(enlist_p1_a, positions_a(2, :), enlist_p2_a, interpolated)
  field_a = field_a * interpolated
  deallocate(interpolated)
  parallel_interpolation_time = mpi_wtime() - t0
  
  
  ! Read the target mesh partition
  t0 = mpi_wtime()
  call read_node(trim(basename_b) // "_" // trim(size_chr) // "_" // trim(rank_chr) // ".node", dim = 2, coords = positions_b)
  call read_ele(trim(basename_b) // "_" // trim(size_chr) // "_" // trim(rank_chr) // ".ele", dim = 2, enlist = enlist_p1_b)
  nnodes_p1_b = size(positions_b, 2)
  nelements_b = size(enlist_p1_b, 2)
  ! Read target mesh halo data ...
  call read_halo(trim(basename_b) // "_" // trim(size_chr), halo_b, level = 2)
  parallel_read_time = parallel_read_time + mpi_wtime() - t0

  t0 = mpi_wtime()
  ! ... and determine the target mesh element ownership
  allocate(ele_owner_b(nelements_b))
  call element_ownership(nnodes_p1_b, enlist_p1_b, halo_b, ele_owner_b)
  ! Generate the donor P2 element-node graph
  allocate(enlist_p2_b(6, nelements_b))
  call p2_connectivity(nnodes_p1_b, enlist_p1_b, nnodes_p2_b, enlist_p2_b)
  ! Construct a target P2 field equal to: x^2
  allocate(field_b(nnodes_p2_b))
  call interpolate_p1_p2(enlist_p1_b, positions_b(1, :), enlist_p2_b, field_b)
  field_b = field_b * field_b
  parallel_interpolation_time = parallel_interpolation_time + mpi_wtime() - t0

  if (rank == 0) print"(a,i15,a,i15)", " A:", nelements_a, ", B:", nelements_b

  t0 = mpi_wtime()
  ! Perform multi-mesh integration
  area_parallel = 0.0D0
  integral_parallel = 0.0D0
  call parallel_supermesh(positions_b, enlist_p1_b, ele_owner_b, &
                        & positions_a, enlist_p1_a, ele_owner_a, &
                        & donor_ele_data, unpack_data_a, intersection_calculation, &
                        & comm = MPI_COMM_WORLD)
  ! Deallocate any remaining unpacked communicated data
  call cleanup_unpack_data_a()
  parallel_time = mpi_wtime() - t0

  ! Sum all process contributions to the multi-mesh integrals
  call MPI_Allreduce(area_parallel, real_buffer, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  area_parallel = real_buffer
  call MPI_Allreduce(integral_parallel, real_buffer, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
  integral_parallel = real_buffer

  ! Sum, Min and Max of parallel read and compute durations
  call MPI_Allreduce(parallel_time, parallel_time_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
  call MPI_Allreduce(parallel_time, parallel_time_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
  call MPI_Allreduce(parallel_time, parallel_time_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
  call MPI_Allreduce(parallel_read_time, parallel_time_read_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
  call MPI_Allreduce(parallel_read_time, parallel_time_read_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
  call MPI_Allreduce(parallel_read_time, parallel_time_read_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
  call MPI_Allreduce(parallel_interpolation_time, parallel_time_interpolation_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
  call MPI_Allreduce(parallel_interpolation_time, parallel_time_interpolation_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
  call MPI_Allreduce(parallel_interpolation_time, parallel_time_interpolation_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr); CHKERRQ(ierr)
#if PROFILE == 1
  call get_times(all_to_all_max, all_to_all_min, all_to_all_sum, &
                 point_to_point_max, point_to_point_min, point_to_point_sum)
#endif

  if(rank == 0) then
    print "(a,f20.15)", "(MIN) Time, parallel = ", parallel_time_min
    print "(a,f20.15)", "(MAX) Time, parallel = ", parallel_time_max
    print "(a,f20.15)", "(SUM) Time, parallel = ", parallel_time_sum
    print "(a,f20.15)", "(AVG) Time, parallel = ", parallel_time_sum / nprocs
    print "(a)" , ""
    print "(a,f20.15)", "(MIN) Interpolation Time, parallel = ", parallel_time_interpolation_min
    print "(a,f20.15)", "(MAX) Interpolation Time, parallel = ", parallel_time_interpolation_max
    print "(a,f20.15)", "(SUM) Interpolation Time, parallel = ", parallel_time_interpolation_sum
    print "(a,f20.15)", "(AVG) Interpolation Time, parallel = ", parallel_time_interpolation_sum / nprocs
    print "(a)" , ""
    print "(a,f20.15)", "(MIN) Read Time, parallel = ", parallel_time_read_min
    print "(a,f20.15)", "(MAX) Read Time, parallel = ", parallel_time_read_max
    print "(a,f20.15)", "(SUM) Read Time, parallel = ", parallel_time_read_sum
    print "(a,f20.15)", "(AVG) Read Time, parallel = ", parallel_time_read_sum / nprocs
    print "(a)" , ""
#if PROFILE == 1
    print "(a,f20.15)", "(MIN) All to all comms = ", all_to_all_min
    print "(a,f20.15)", "(MAX) All to all comms = ", all_to_all_max
    print "(a,f20.15)", "(SUM) All to all comms = ", all_to_all_sum
    print "(a,f20.15)", "(AVG) All to all comms = ", all_to_all_sum / nprocs
    print "(a,f20.15)", "(MIN) Point to point comms = ", point_to_point_min
    print "(a,f20.15)", "(MAX) Point to point comms = ", point_to_point_max
    print "(a,f20.15)", "(SUM) Point to point comms = ", point_to_point_sum
    print "(a,f20.15)", "(AVG) Point to point comms = ", point_to_point_sum / nprocs
    print "(a)" , ""
#endif
    ! Display the multi-mesh integrals on rank 0
    print "(a,e26.18e3,a,e26.18e3,a)", "Area     = ", area_parallel, " (error = ", abs(area_parallel - 0.5D0), ")"
    print "(a,e26.18e3,a,e26.18e3,a)", "Integral = ", integral_parallel, " (error = ", abs(integral_parallel - 2.7083333333333272D-02), ")"
  end if

  ! Cleanup
  deallocate(positions_a, enlist_p1_a, ele_owner_a, enlist_p2_a, field_a, &
           & positions_b, enlist_p1_b, ele_owner_b, enlist_p2_b, field_b)
  call deallocate(halo_a)
  call deallocate(halo_b)

contains

  ! Generate a P2 element-node graph
  subroutine p2_connectivity(nnodes_p1, enlist_p1, nnodes_p2, enlist_p2)
    ! Number of P1 nodes
    integer, intent(in) :: nnodes_p1
    ! P1 element-node graph
    ! Shape: 3 x nelements
    integer, dimension(:, :), intent(in) :: enlist_p1
    ! Number of P2 nodes
    integer, intent(out) :: nnodes_p2
    ! P2 element-node graph
    ! Shape: 6 x nelements
    integer, dimension(:, :), intent(out) :: enlist_p2
    
    integer :: ele, lnode, nelements, node_p1, node_p1_1, node_p1_2, node_p2
    type(integer_hash_table) :: node_map_p1_p2
    type(integer_hash_table), dimension(:), allocatable :: node_map_p2
    
    nelements = size(enlist_p1, 2)
    
    call allocate(node_map_p1_p2)
    allocate(node_map_p2(nnodes_p1))
    do node_p1 = 1, nnodes_p1
      call allocate(node_map_p2(node_p1))
    end do
    nnodes_p2 = 0
        
    do ele = 1, nelements    
      do lnode = 1, 6
        select case(lnode)
          case(1:3)
            node_p1 = enlist_p1(lnode, ele)
            if(has_key(node_map_p1_p2, node_p1)) then
              node_p2 = fetch(node_map_p1_p2, node_p1)
            else
              node_p2 = nnodes_p2 + 1
              nnodes_p2 = nnodes_p2 + 1
              call insert(node_map_p1_p2, node_p1, node_p2)
            end if
          case(4)
            node_p1_1 = enlist_p1(1, ele)
            node_p1_2 = enlist_p1(2, ele)
            node_p1 = min(node_p1_1, node_p1_2)
            node_p1_2 = max(node_p1_1, node_p1_2)
            node_p1_1 = node_p1
            if(has_key(node_map_p2(node_p1_1), node_p1_2)) then
              node_p2 = fetch(node_map_p2(node_p1_1), node_p1_2)
            else
              node_p2 = nnodes_p2 + 1
              nnodes_p2 = nnodes_p2 + 1
              call insert(node_map_p2(node_p1_1), node_p1_2, node_p2)
            end if
          case(5)
            node_p1_1 = enlist_p1(2, ele)
            node_p1_2 = enlist_p1(3, ele)
            node_p1 = min(node_p1_1, node_p1_2)
            node_p1_2 = max(node_p1_1, node_p1_2)
            node_p1_1 = node_p1
            if(has_key(node_map_p2(node_p1_1), node_p1_2)) then
              node_p2 = fetch(node_map_p2(node_p1_1), node_p1_2)
            else
              node_p2 = nnodes_p2 + 1
              nnodes_p2 = nnodes_p2 + 1
              call insert(node_map_p2(node_p1_1), node_p1_2, node_p2)
            end if
          case(6)
            node_p1_1 = enlist_p1(1, ele)
            node_p1_2 = enlist_p1(3, ele)
            node_p1 = min(node_p1_1, node_p1_2)
            node_p1_2 = max(node_p1_1, node_p1_2)
            node_p1_1 = node_p1
            if(has_key(node_map_p2(node_p1_1), node_p1_2)) then
              node_p2 = fetch(node_map_p2(node_p1_1), node_p1_2)
            else
              node_p2 = nnodes_p2 + 1
              nnodes_p2 = nnodes_p2 + 1
              call insert(node_map_p2(node_p1_1), node_p1_2, node_p2)
            end if
          case default
            stop 1
        end select
        enlist_p2(lnode, ele) = node_p2
      end do      
    end do
    
    call deallocate(node_map_p1_p2)
    do node_p1 = 1, nnodes_p1
      call deallocate(node_map_p2(node_p1))
    end do
    deallocate(node_map_p2)
    
  end subroutine p2_connectivity
  
  ! Interpolate a P1 field onto a P2 function space
  subroutine interpolate_p1_p2(enlist_p1, field_p1, enlist_p2, field_p2)
    ! P1 element-node graph
    ! Shape: 3 x nelements
    integer, dimension(:, :), intent(in) :: enlist_p1
    ! P1 field
    ! Shape: nnodes_p1
    real, dimension(:), intent(in) :: field_p1
    ! P2 element-node graph
    ! Shape: 6 x nelements
    integer, dimension(:, :), intent(in) :: enlist_p2
    ! P2 field
    ! Shape: nnodes_p2
    real, dimension(:), intent(out) :: field_p2
    
    integer :: ele, lnode, nelements, node_p2, nnodes_p2
    logical, dimension(:), allocatable :: seen_node_p2
    
    nelements = size(enlist_p2, 2)
    nnodes_p2 = size(field_p2)
    
    allocate(seen_node_p2(nnodes_p2))
    seen_node_p2 = .false.
    
    do ele = 1, nelements
      do lnode = 1, 6
        node_p2 = enlist_p2(lnode, ele)
        if(seen_node_p2(node_p2)) cycle
        seen_node_p2(node_p2) = .true.
        select case(lnode)
          case(1:3)
            field_p2(node_p2) = field_p1(enlist_p1(lnode, ele))
          case(4)
            field_p2(node_p2) = 0.5D0 * (field_p1(enlist_p1(1, ele)) + field_p1(enlist_p1(2, ele)))
          case(5)
            field_p2(node_p2) = 0.5D0 * (field_p1(enlist_p1(2, ele)) + field_p1(enlist_p1(3, ele)))
          case(6)
            field_p2(node_p2) = 0.5D0 * (field_p1(enlist_p1(1, ele)) + field_p1(enlist_p1(3, ele)))
          case default
            stop 1
        end select
      end do
    end do
    
    deallocate(seen_node_p2)
  
  end subroutine interpolate_p1_p2

  ! Given the provided mesh vertices and elements, pack data for communication
  subroutine donor_ele_data(nodes_a, eles_a, data_a)
    ! Mesh vertices to be communicated
    integer, dimension(:), intent(in) :: nodes_a
    ! Mesh elements to be communicated
    integer, dimension(:), intent(in) :: eles_a
    ! Packed data for communication
    integer(kind = c_int8_t), dimension(:), allocatable, intent(out) :: data_a
 
    integer :: ele, i, lnode, ndata_a, node_a, nnodes_p2_a, position
    integer, dimension(:, :), allocatable :: data_enlist_p2_a
    type(integer_hash_table) :: node_map
    type(integer_set) :: nodes_p2_a
    
    ! For which P2 nodes do we need to send data?
    call allocate(nodes_p2_a)
    do i = 1, size(eles_a)
      ele = eles_a(i)
      do lnode = 1, 6
        node_a = enlist_p2_a(lnode, ele)
        call insert(nodes_p2_a, node_a)
      end do
    end do
    nnodes_p2_a = key_count(nodes_p2_a)
    
    ! Construct a map from communicated nodes to local nodes ...
    call allocate(node_map)
    do i = 1, nnodes_p2_a
      node_a = fetch(nodes_p2_a, i)
      call insert(node_map, node_a, i)
    end do
    ! ... and use this to construct the communicated P2 element-node graph
    allocate(data_enlist_p2_a(6, size(eles_a)))
    do i = 1, size(eles_a)
      ele = eles_a(i)
      do lnode = 1, 6
        node_a = enlist_p2_a(lnode, ele)
        data_enlist_p2_a(lnode, i) = fetch(node_map, node_a)
      end do
    end do
    
    ! Gather P2 field values for communication
    allocate(data_field_a(nnodes_p2_a))
    do i = 1, nnodes_p2_a
      node_a = fetch(nodes_p2_a, i)
      data_field_a(i) = field_a(node_a)
    end do
    
    ! Pack data for communication:
    !   1 integer                         -- number of elements
    !   1 integer                         -- number of P2 nodes
    !   (6 x number of elements) integers -- communicated P2 element-node graph
    !   (number of P2 nodes) reals        -- communicated P2 field values
    ndata_a = (2 + 6 * size(eles_a)) * integer_extent + nnodes_p2_a * real_extent
    allocate(data_a(ndata_a))
    position = 0
    call MPI_Pack(size(eles_a), 1, MPI_INTEGER, data_a, ndata_a, position, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    call MPI_Pack(nnodes_p2_a, 1, MPI_INTEGER, data_a, ndata_a, position, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    call MPI_Pack(data_enlist_p2_a, 6 * size(eles_a), MPI_INTEGER, data_a, ndata_a, position, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    call MPI_Pack(data_field_a, nnodes_p2_a, MPI_DOUBLE_PRECISION, data_a, ndata_a, position, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    
    call deallocate(nodes_p2_a)
    call deallocate(node_map)
    deallocate(data_enlist_p2_a, data_field_a)
    
  end subroutine donor_ele_data
  
  ! Unpack communicated data
  subroutine unpack_data_a(data_a)
    ! Packed communicated data
    integer(kind = c_int8_t), dimension(:), intent(in) :: data_a
    
    integer :: position
    
    ! Deallocate any previously unpacked communicated data
    call cleanup_unpack_data_a()
    
    position = 0
    ! Unpack the number of elements
    call MPI_Unpack(data_a, size(data_a), position, data_nelements_a, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    ! Unpack the number of P2 nodes
    call MPI_Unpack(data_a, size(data_a), position, data_nnodes_p2_a, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    ! Unpack the P2 element-node graph
    allocate(data_enlist_p2_a(6, data_nelements_a))
    call MPI_Unpack(data_a, size(data_a), position, data_enlist_p2_a, 6 * data_nelements_a, MPI_INTEGER, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    ! Unpack the P2 field values
    allocate(data_field_a(data_nnodes_p2_a))
    call MPI_Unpack(data_a, size(data_a), position, data_field_a, data_nnodes_p2_a, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr);  CHKERRQ(ierr)
    
  end subroutine unpack_data_a
  
  ! Deallocate any previously unpacked communicated data
  subroutine cleanup_unpack_data_a()
    if(allocated(data_enlist_p2_a)) deallocate(data_enlist_p2_a)
    if(allocated(data_field_a)) deallocate(data_field_a)
  end subroutine cleanup_unpack_data_a

  ! Evaluate P2 basis functions at a given point
  pure function basis_functions_p2(cell_coords, coord) result(fns)
    ! Triangle vertex coordinates
    ! Shape: dim x loc_p1
    real, dimension(2, 3), intent(in) :: cell_coords
    ! Coordinate at which to evaluate the basis functions
    ! Shape: dim
    real, dimension(2), intent(in) :: coord

    ! Basis function values
    ! Shape: loc_p2
    real, dimension(6) :: fns

    real :: x, y
    real, dimension(2) :: e1, e2, lcoord
    real, dimension(2, 2) :: jac

    e1 = cell_coords(:, 2) - cell_coords(:, 1)
    e2 = cell_coords(:, 3) - cell_coords(:, 1)

    jac(1, 1) =  e2(2);  jac(1, 2) = -e2(1)
    jac(2, 1) = -e1(2);  jac(2, 2) =  e1(1)
    jac = jac / (jac(1, 1) * jac(2, 2) - jac(1, 2) * jac(2, 1))

    lcoord = matmul(jac, coord - cell_coords(:, 1))
    x = lcoord(1)
    y = lcoord(2)
    
    fns(1) = (1.0D0 - x - y) * (1.0D0 - 2.0D0 * x - 2.0D0 * y)
    fns(2) = x * (2.0D0 * x - 1.0D0)
    fns(3) = y * (2.0D0 * y - 1.0D0)
    fns(4) = 4.0D0 * x * (1.0D0 - x - y)
    fns(5) = 4.0D0 * x * y
    fns(6) = 4.0D0 * y * (1.0D0 - x - y)

  end function basis_functions_p2

  ! Interpolate a P2 function at given point
  pure function interpolate_p2(cell_coords_d, cell_x_d, coord_s) result(x_s)
    ! Triangle vertex coordinates
    ! Shape: dim x loc_p1
    real, dimension(2, 3), intent(in) :: cell_coords_d
    ! P2 nodal values
    ! Shape: loc_p2
    real, dimension(6), intent(in) :: cell_x_d
    ! Coordinate at which to evaluate the P2 function
    ! Shape: dim
    real, dimension(2), intent(in) :: coord_s

    real :: x_s

    x_s = dot_product(basis_functions_p2(cell_coords_d, coord_s), cell_x_d)

  end function interpolate_p2
  
  ! Perform calculations on the local supermesh
  subroutine intersection_calculation(positions_b, positions_a, positions_c, nodes_a, ele_b, ele_a, local)
    ! Target mesh element vertex coordinates
    ! Shape: dim x loc_b
    real, dimension(:, :), intent(in) :: positions_b
    ! Donor mesh element vertex coordinates
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
    
    integer :: ele_c, lnode
    integer, dimension(:, :), pointer :: lenlist_p2_a
    real :: area
    real, dimension(6) :: field_a_c, field_b_c
    real, dimension(:), pointer :: lfield_a
    
    if(local) then
      ! If this is a local calculation, use the local P2 element-node graph and
      ! field data
      lenlist_p2_a => enlist_p2_a
      lfield_a => field_a
    else
      ! Otherwise, use the unpacked communicated element-node graph and P2 field
      ! data
      lenlist_p2_a => data_enlist_p2_a
      lfield_a => data_field_a
    end if
  
    do ele_c = 1, size(positions_c, 3)
      ! Compute the supermesh triangle area
      area = triangle_area(positions_c(:, :, ele_c))
      ! Local contribution to the intersection area
      area_parallel = area_parallel + area
      ! Interpolate the donor and target P2 functions onto a P2 space on the
      ! supermesh element
      do lnode = 1, 6
        select case(lnode)
          case(1:3)
            field_a_c(lnode) = interpolate_p2(positions_a, lfield_a(lenlist_p2_a(:, ele_a)), positions_c(:, lnode, ele_c))
            field_b_c(lnode) = interpolate_p2(positions_b, field_b(enlist_p2_b(:, ele_b)), positions_c(:, lnode, ele_c))
          case(4)
            field_a_c(lnode) = interpolate_p2(positions_a, lfield_a(lenlist_p2_a(:, ele_a)), 0.5D0 * (positions_c(:, 1, ele_c) + positions_c(:, 2, ele_c)))
            field_b_c(lnode) = interpolate_p2(positions_b, field_b(enlist_p2_b(:, ele_b)), 0.5D0 * (positions_c(:, 1, ele_c) + positions_c(:, 2, ele_c)))
          case(5)
            field_a_c(lnode) = interpolate_p2(positions_a, lfield_a(lenlist_p2_a(:, ele_a)), 0.5D0 * (positions_c(:, 2, ele_c) + positions_c(:, 3, ele_c)))
            field_b_c(lnode) = interpolate_p2(positions_b, field_b(enlist_p2_b(:, ele_b)), 0.5D0 * (positions_c(:, 2, ele_c) + positions_c(:, 3, ele_c)))
          case(6)
            field_a_c(lnode) = interpolate_p2(positions_a, lfield_a(lenlist_p2_a(:, ele_a)), 0.5D0 * (positions_c(:, 1, ele_c) + positions_c(:, 3, ele_c)))
            field_b_c(lnode) = interpolate_p2(positions_b, field_b(enlist_p2_b(:, ele_b)), 0.5D0 * (positions_c(:, 1, ele_c) + positions_c(:, 3, ele_c)))
          case default
            stop 1
        end select
      end do
      ! Local contribution to the multi-mesh inner product
      integral_parallel = integral_parallel + 2.0D0 * area * dot_product(field_a_c, matmul(mass_p2, field_b_c))
    end do
    
  end subroutine intersection_calculation
  
end subroutine benchmark_parallel_p2_inner_product
