subroutine benchmark_parallel

  use libsupermesh_unittest_tools
  use libsupermesh_construction
  use libsupermesh_fields
  use libsupermesh_read_triangle
  use libsupermesh_tri_intersection_module
  use libsupermesh_unittest_tools
  use libsupermesh_intersection_finder
  
  implicit none

#include <finclude/petsc.h90>

  integer :: i, nnodes, ele_A, ele_B, ele_C, n_trisC, mpi_num_procs, mpi_my_id, mpi_my_error
  character(len=5) :: mpi_my_id_character
  type(tri_type) :: tri_A, tri_B
  type(tri_type), dimension(tri_buf_size) :: trisC
  logical :: pass
  real :: t_0, t1, t2, area_serial = 0.0, area_parallel = 0.0, & 
       & time_serial = 0.0, time_parallel = 0.0, mesh_serial = 0.0, mesh_parallel = 0.0
  real, dimension(:), allocatable :: areas_parallel, times_parallel
  type(vector_field) :: positionsA, positionsB
  integer, parameter :: dim = 2, loc = 3
  character(len=9999) :: filenameA, filenameB
  integer, parameter :: mpi_my_root = 0
  integer, dimension(:), allocatable :: nodes
  real, parameter :: tol = 1.0e3 * epsilon(0.0)
  logical :: fail = .FALSE.
  real, dimension(:,:), allocatable :: bbox_a, bbox_b
  PetscErrorCode :: ierr
  
!  call MPI_INIT ( mpi_my_error ); CHKERRQ(mpi_my_error)

  !     find out MY process ID, and how many processes were started.
  call MPI_COMM_RANK (MPI_COMM_WORLD, mpi_my_id, mpi_my_error); CHKERRQ(mpi_my_error)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, mpi_num_procs, mpi_my_error); CHKERRQ(mpi_my_error)

  call PETScInitialize(PETSC_NULL_CHARACTER, ierr);  CHKERRQ(ierr)
  pass = .true.

! Do everything in serial  
  if ( mpi_my_id .eq. 0 ) then  
    t_0 = mpi_wtime()
    positionsA = read_triangle_files("data/square_0_2", dim)
    positionsB = read_triangle_files("data/square_0_1", dim)
    mesh_serial = mpi_wtime() - t_0

    print "(a,e25.17e3)", "Mesh input time  = ", mesh_serial
    print "(a,i10,a,i10,a)", "Element Count (tris) A:",ele_count(positionsA),", Node Count (cells) A:",node_count(positionsA),"."
    print "(a,i10,a,i10,a)", "Element Count (tris) B:",ele_count(positionsB),", Node Count (cells) B:",node_count(positionsB),"."

    t1 = mpi_wtime()
    do ele_A=1,ele_count(positionsA)
      tri_A%v = ele_val(positionsA, ele_A)

      do ele_B=1,ele_count(positionsB) 
        ! B. Use the libSuperMesh internal triangle intersector (using derived types as input)
        tri_B%v = ele_val(positionsB, ele_B)

        call intersect_tris(tri_A, tri_B, trisC, n_trisC)

        do ele_C=1,n_trisC
          area_serial = area_serial + triangle_area(trisC(ele_C)%v)
        end do
      end do
    end do
    t2 = mpi_wtime()
    time_serial = t2 - t1

    print "(a,e25.17e3)", "Serial intersection time  = ", time_serial
    print "(a,e25.17e3,a)", "Serial intersection area:", area_serial,"."
    
    call deallocate(positionsA)
    call deallocate(positionsB)
  end if
  
!  ToDo TODO todo FIX HACK
!  Remove the comments when running non debug code
!  if ( mpi_my_id .eq. 0 ) then 
    allocate(areas_parallel(0:mpi_num_procs-1), times_parallel(0:mpi_num_procs-1))
!  end if
!  ToDo TODO todo FIX HACK
!  Remove the next lines
  call MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
  print "(a)", "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
  call FLUSH()
  call MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
  
!  ToDo TODO todo FIX HACK
!  Refactor the next lines. We only load this file here
!  because it messes up our printouts.
!  do i=0,mpi_num_procs
!    if ( mpi_my_id .eq. i ) then 
!      positionsA = read_triangle_files("data/square_0_2", dim)
!!      positionsB = read_triangle_files("data/square_0_1", dim)
!      call FLUSH()
!    end if
!  end do
!  ToDo TODO todo FIX HACK
!  Remove the next lines
!  call MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
!  print "(a)", "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
!  call FLUSH()
!  call MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)

!  ToDo TODO todo FIX HACK
!  Remove the DO loop and the IF check
  do i=0,mpi_num_procs
    call FLUSH()
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
    if ( mpi_my_id .eq. i ) then 
      area_parallel = 0.0
      write(mpi_my_id_character,'(I5)') mpi_my_id

      t_0 = mpi_wtime()
      filenameA = trim(adjustl("data/square_0_2_"))//trim(adjustl(mpi_my_id_character))

      call read_halo_files(trim(filenameA), nnodes, nodes)
!  ToDo TODO todo FIX HACK
      print "(a,i5,a,i10,a,a,a)", "MPI Process:",mpi_my_id,", has nodes:",nnodes,", read from:",trim(filenameA),"."
      call FLUSH()
      if (nnodes < dim + 1) cycle

!      positionsA = read_triangle_files("data/square_0_2", dim)
      positionsA = read_triangle_files(trim(filenameA), dim, nnodes = nnodes, nodes = nodes)



      nnodes = 0
      deallocate(nodes)
      filenameB = trim(adjustl("data/square_0_1_"))//trim(adjustl(mpi_my_id_character))

!      call read_halo_files(trim(filenameB), nnodes, nodes)
!  ToDo TODO todo FIX HACK
!      print "(a,i5,a,i10,a,a,a)", "MPI Process:",mpi_my_id,", has nodes:",nnodes,", read from:",trim(filenameB),"."
      call FLUSH()
!      if (nnodes < dim + 1) cycle

      positionsB = read_triangle_files("data/square_0_1", dim)
!      positionsB = read_triangle_files(trim(filenameB), dim, nnodes = nnodes, nodes = nodes)

      mesh_parallel = mpi_wtime() - t_0

!  ToDo TODO todo FIX HACK
      print "(i5,a,e25.17e3)", mpi_my_id,": Mesh input time  = ", mesh_parallel
      print "(i5,a,i10,a,i10,a)", mpi_my_id,": Element Count (tris) A:",ele_count(positionsA),", Node Count (cells) A:",node_count(positionsA),"."
      print "(i5,a,i10,a,i10,a)", mpi_my_id,": Element Count (tris) B:",ele_count(positionsB),", Node Count (cells) B:",node_count(positionsB),"."

      allocate(bbox_a(positionsA%dim,2))
      bbox_a = bbbox(positionsA)
      allocate(bbox_b(positionsB%dim,2))
      bbox_b = bbbox(positionsB)
      write(*,*) "bbox_a                                :",bbox_a,"."
      write(*,*) "bbox_b                                :",bbox_b,"."
!  ToDo TODO todo FIX HACK
      print "(a)", "======================================================================"
      call FLUSH()

      t1 = mpi_wtime()
      do ele_A=1,ele_count(positionsA)
!        write(*,*) "ele_A:",ele_A,"."
!        write(*,*) "positionsA%mesh%ndglno(ele_A) :",positionsA%mesh%ndglno(ele_A),"."
!        write(*,*) "ele_val(",ele_A,") :",ele_val(positionsA, ele_A),"."
!        write(*,*) "positionsA%val(:, ele_A) :",positionsA%val(:, :),"."
        call FLUSH()
        tri_A%v = ele_val(positionsA, ele_A)
        do ele_B=1,ele_count(positionsB) 
          ! B. Use the libSuperMesh internal triangle intersector (using derived types as input)
          tri_B%v = ele_val(positionsB, ele_B)

          call intersect_tris(tri_A, tri_B, trisC, n_trisC)

          do ele_C=1,n_trisC
            area_parallel = area_parallel + triangle_area(trisC(ele_C)%v)
          end do
        end do
      end do
      t2 = mpi_wtime()
      time_parallel = t2 - t1

      print "(i5,a,e25.17e3)", mpi_my_id, ": intersection time  = ", time_parallel
      print "(i5,a,e25.17e3,a)", mpi_my_id, ": intersection area:", area_parallel,"."
    
      call deallocate(positionsA)
      call deallocate(positionsB)
!  ToDo TODO todo FIX HACK
      print "(a)", "======================================================================"
      print "(a)", "======================================================================"
    end if
    call FLUSH()
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
  end do
  
  call FLUSH()
  call MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
  call FLUSH()

  ! Gather remote results:
  call MPI_Gather(area_parallel, 1, MPI_DOUBLE_PRECISION, &
      & areas_parallel, 1, MPI_DOUBLE_PRECISION, &
      & mpi_my_root, MPI_COMM_WORLD, mpi_my_error)
  call MPI_Gather(time_parallel, 1, MPI_DOUBLE_PRECISION, &
      & times_parallel, 1, MPI_DOUBLE_PRECISION, &
      & mpi_my_root, MPI_COMM_WORLD, mpi_my_error)

  if ( mpi_my_id .eq. 0 ) then 
    write (*,*) times_parallel
    write (*,*) areas_parallel
    print "(i5,a,e25.17e3,a)", mpi_my_id, ": Total parallel intersection time:", sum(times_parallel),"."
    print "(i5,a,e25.17e3,a)", mpi_my_id, ": Total parallel intersection area:", sum(areas_parallel),"."
    fail = fnequals(sum(areas_parallel), area_serial, tol = tol)
    call report_test("[benchmark_parallel areas]", fail, .false., "Should give the same areas of intersection")
  end if
  
  deallocate(areas_parallel, times_parallel, bbox_a, bbox_b)
  
  
  
  
  
  


  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
!  ToDo TODO todo FIX HACK
!  Remove the comments when running non debug code
!  if ( mpi_my_id .eq. 0 ) then 
    allocate(areas_parallel(0:mpi_num_procs-1), times_parallel(0:mpi_num_procs-1))
!  end if
!  ToDo TODO todo FIX HACK
!  Remove the next lines
  call MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
  print "(a)", "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
  call FLUSH()
  call MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
  
!  ToDo TODO todo FIX HACK
!  Remove the DO loop and the IF check
  do i=0,mpi_num_procs
    call FLUSH()
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
    if ( mpi_my_id .eq. i ) then 
      area_parallel = 0.0
      write(mpi_my_id_character,'(I5)') mpi_my_id

      t_0 = mpi_wtime()
      filenameA = trim(adjustl("data/square_0_2_"))//trim(adjustl(mpi_my_id_character))

      call read_halo_files(trim(filenameA), nnodes, nodes)
!  ToDo TODO todo FIX HACK
      print "(a,i5,a,i10,a,a,a)", "MPI Process:",mpi_my_id,", has nodes:",nnodes,", read from:",trim(filenameA),"."
      call FLUSH()
      if (nnodes < dim + 1) cycle

!      positionsA = read_triangle_files("data/square_0_2", dim)
      positionsA = read_triangle_files(trim(filenameA), dim, nnodes = nnodes, nodes = nodes)



      nnodes = 0
      deallocate(nodes)
      filenameB = trim(adjustl("data/square_0_1_"))//trim(adjustl(mpi_my_id_character))

      call read_halo_files(trim(filenameB), nnodes, nodes)
!  ToDo TODO todo FIX HACK
      print "(a,i5,a,i10,a,a,a)", "MPI Process:",mpi_my_id,", has nodes:",nnodes,", read from:",trim(filenameB),"."
      call FLUSH()
      if (nnodes < dim + 1) cycle

!      positionsB = read_triangle_files("data/square_0_1", dim)
      positionsB = read_triangle_files(trim(filenameB), dim, nnodes = nnodes, nodes = nodes)

      mesh_parallel = mpi_wtime() - t_0

!  ToDo TODO todo FIX HACK
      print "(i5,a,e25.17e3)", mpi_my_id,": Mesh input time  = ", mesh_parallel
      print "(i5,a,i10,a,i10,a)", mpi_my_id,": Element Count (tris) A:",ele_count(positionsA),", Node Count (cells) A:",node_count(positionsA),"."
      print "(i5,a,i10,a,i10,a)", mpi_my_id,": Element Count (tris) B:",ele_count(positionsB),", Node Count (cells) B:",node_count(positionsB),"."
      
      allocate(bbox_a(positionsA%dim,2))
      bbox_a = bbbox(positionsA)
      allocate(bbox_b(positionsB%dim,2))
      bbox_b = bbbox(positionsB)
      write(*,*) "bbox_a                                :",bbox_a,"."
      write(*,*) "bbox_b                                :",bbox_b,"."
!  ToDo TODO todo FIX HACK
      print "(a)", "======================================================================"
      call FLUSH()

    end if
  end do

  call FLUSH()
  call MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
  call FLUSH()
  
  ! Distribute partitionA to all nodes
!  call MPI_Gather(area_parallel, 1, MPI_DOUBLE_PRECISION, &
!      & areas_parallel, 1, MPI_DOUBLE_PRECISION, &
!      & mpi_my_root, MPI_COMM_WORLD, mpi_my_error)

  call FLUSH()
  call MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
  call FLUSH()

  ! Compute the intersection
  do i=0,mpi_num_procs
    call FLUSH()
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
    if ( mpi_my_id .eq. i ) then 
      t1 = mpi_wtime()

      do ele_A=1,ele_count(positionsA)
!        write(*,*) "ele_A:",ele_A,"."
!        write(*,*) "positionsA%mesh%ndglno(ele_A) :",positionsA%mesh%ndglno(ele_A),"."
        write(*,*) "ele_val(",ele_A,") :",ele_val(positionsA, ele_A),"."
!        write(*,*) "positionsA%val(:, ele_A) :",positionsA%val(:, :),"."
        call FLUSH()
        tri_A%v = ele_val(positionsA, ele_A)
        do ele_B=1,ele_count(positionsB) 
          ! B. Use the libSuperMesh internal triangle intersector (using derived types as input)
          tri_B%v = ele_val(positionsB, ele_B)

          call intersect_tris(tri_A, tri_B, trisC, n_trisC)

          do ele_C=1,n_trisC
            area_parallel = area_parallel + triangle_area(trisC(ele_C)%v)
          end do
        end do
      end do
      t2 = mpi_wtime()
      time_parallel = t2 - t1

      print "(i5,a,e25.17e3)", mpi_my_id, ": intersection time  = ", time_parallel
      print "(i5,a,e25.17e3,a)", mpi_my_id, ": intersection area:", area_parallel,"."
    
      call deallocate(positionsA)
      call deallocate(positionsB)
!  ToDo TODO todo FIX HACK
      print "(a)", "======================================================================"
      print "(a)", "======================================================================"
    end if
    call FLUSH()
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
  end do
  
  call FLUSH()
  call MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
  call FLUSH()

  ! Gather remote results:
  call MPI_Gather(area_parallel, 1, MPI_DOUBLE_PRECISION, &
      & areas_parallel, 1, MPI_DOUBLE_PRECISION, &
      & mpi_my_root, MPI_COMM_WORLD, mpi_my_error)
  call MPI_Gather(time_parallel, 1, MPI_DOUBLE_PRECISION, &
      & times_parallel, 1, MPI_DOUBLE_PRECISION, &
      & mpi_my_root, MPI_COMM_WORLD, mpi_my_error)

  if ( mpi_my_id .eq. 0 ) then 
    write (*,*) times_parallel
    write (*,*) areas_parallel
    print "(i5,a,e25.17e3,a)", mpi_my_id, ": Total parallel intersection time:", sum(times_parallel),"."
    print "(i5,a,e25.17e3,a)", mpi_my_id, ": Total parallel intersection area:", sum(areas_parallel),"."
    fail = fnequals(sum(areas_parallel), area_serial, tol = tol)
    call report_test("[benchmark_parallel areas]", fail, .false., "Should give the same areas of intersection")
  end if
  
  deallocate(areas_parallel, times_parallel, bbox_a, bbox_b)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

  call PETScFinalize(ierr); CHKERRQ(mpi_my_error)
!  call MPI_FINALIZE ( mpi_my_error ); CHKERRQ(mpi_my_error)
  if(.not. pass) stop 1

end subroutine benchmark_parallel