subroutine benchmark_parallel

  use libsupermesh_unittest_tools
  use libsupermesh_construction
  use libsupermesh_fields_dummy
  use libsupermesh_read_triangle_2
  use libsupermesh_tri_intersection_module
  use libsupermesh_unittest_tools
  
  implicit none

#include <finclude/petsc.h90>

  integer :: i, nnodes, ele_A, ele_B, ele_C, n_trisC, mpi_num_procs, mpi_my_id, mpi_my_error
  character(len=5) :: mpi_my_id_character
  type(tri_type) :: tri_A, tri_B
  type(tri_type), dimension(tri_buf_size) :: trisC
  logical :: pass
  real :: t_0, t1, t2, area_serial = 0.0, area_parallel = 0.0
  real, dimension(:), allocatable :: areas_parallel
  type(vector_field) :: positionsA, positionsB
  integer, parameter :: dim = 2, loc = 3
  character(len=9999) :: filename
  integer, parameter :: mpi_my_root = 0
  PetscErrorCode :: ierr
  
  call MPI_INIT ( mpi_my_error )

  !     find out MY process ID, and how many processes were started.
  call MPI_COMM_RANK (MPI_COMM_WORLD, mpi_my_id, mpi_my_error)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, mpi_num_procs, mpi_my_error)

  call PETScInitialize(PETSC_NULL_CHARACTER, ierr);  CHKERRQ(ierr)
  pass = .true.

! Do everything in serial  
  if ( mpi_my_id .eq. 0 ) then  
    t_0 = mpi_wtime()
    positionsA = read_triangle_files("data/square_0_2", dim)
    positionsB = read_triangle_files("data/square_0_1", dim)
    print "(a,e25.17e3)", "Mesh input time  = ", mpi_wtime() - t_0
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

    print "(a,e25.17e3)", "Serial intersection time  = ", t2 - t1
    print "(a,e25.17e3,a)", "Serial intersection area:", area_serial,"."
    
    call deallocate(positionsA)
    call deallocate(positionsB)
  end if
  
!  ToDo TODO todo FIX HACK
!  Remove the comments when running non debug code
!  if ( mpi_my_id .eq. 0 ) then 
    allocate(areas_parallel(0:mpi_num_procs-1))
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
  do i=0,mpi_num_procs
    if ( mpi_my_id .eq. i ) then 
      positionsB = read_triangle_files("data/square_0_1", dim)
      call FLUSH()
    end if
  end do
!  ToDo TODO todo FIX HACK
!  Remove the next lines
  call MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
  print "(a)", "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
  call FLUSH()
  call MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)

!  ToDo TODO todo FIX HACK
!  Remove the DO loop and the IF check
  do i=0,mpi_num_procs
    if ( mpi_my_id .eq. i ) then 
      write(mpi_my_id_character,'(I5)') mpi_my_id
      filename = trim(adjustl("data/square_0_2_"))//trim(adjustl(mpi_my_id_character))
      nnodes = read_halo_files(trim(filename))
!  ToDo TODO todo FIX HACK
      print "(a,i5,a,i10,a,a,a)", "MPI Process:",mpi_my_id,", has nodes:",nnodes,", read from:",trim(filename),"."
  
      t_0 = mpi_wtime()
      positionsA = read_triangle_files(trim(filename), dim, nodes = nnodes)
!  positionsB = read_triangle_files("data/square_0_1", dim)
!  ToDo TODO todo FIX HACK
      print "(i5,a,e25.17e3)", mpi_my_id,": Mesh input time  = ", mpi_wtime() - t_0
      print "(i5,a,i10,a,i10,a)", mpi_my_id,": Element Count (tris) A:",ele_count(positionsA),", Node Count (cells)A:",node_count(positionsA),"."
      print "(i5,a,i10,a,i10,a)", mpi_my_id,": Element Count (tris) B:",ele_count(positionsB),", Node Count (cells) B:",node_count(positionsB),"."

      area_parallel = 0.0
      t1 = mpi_wtime()
      do ele_A=1,ele_count(positionsA)
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

      print "(i5,a,e25.17e3)", mpi_my_id, ": intersection time  = ", t2 - t1
      print "(i5,a,e25.17e3,a)", mpi_my_id, ": intersection area:", area_parallel,"."
    
      call deallocate(positionsA)
      call deallocate(positionsB)
!  ToDo TODO todo FIX HACK
      print "(a)", "======================================================================"
      print "(a)", "======================================================================"
    end if
!  ToDo TODO todo FIX HACK
!  Remove the next lines
    call FLUSH()
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
  end do

  ! Gather remote results:
  call MPI_Gather(area_parallel, 1, MPI_DOUBLE_PRECISION, &
      & areas_parallel, 1, MPI_DOUBLE_PRECISION, &
      & mpi_my_root, MPI_COMM_WORLD, mpi_my_error)

  if ( mpi_my_id .eq. 0 ) then 
    write (*,*) areas_parallel
    print "(i5,a,e25.17e3,a)", mpi_my_id, ": Total intersection area:", sum(areas_parallel),"."
  end if
  
  call PETScFinalize(ierr)
  call MPI_FINALIZE ( mpi_my_error )
  if(.not. pass) stop 1

end subroutine benchmark_parallel