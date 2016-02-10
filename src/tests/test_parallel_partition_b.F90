#include "libsupermesh_debug.h"

subroutine test_parallel_partition_b

  use libsupermesh_unittest_tools
  use libsupermesh_supermesh
  use libsupermesh_fields
  use libsupermesh_read_triangle
  use libsupermesh_tri_intersection
  use libsupermesh_unittest_tools
  use libsupermesh_intersection_finder
  use libsupermesh_read_halos
  use libsupermesh_halo_ownership

  implicit none

#include <finclude/petsc.h90>

  integer :: i, ele_A, ele_B, ele_C, n_trisC, mpi_num_procs, mpi_my_id, mpi_my_error, &
       & serial_ele_A, serial_ele_B, parallel_ele_A, parallel_ele_B, local_sum
  character(len=5) :: mpi_my_id_character, mpi_num_procs_character
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
  real, parameter :: tol = 1.0e3 * epsilon(0.0)
  logical :: fail = .FALSE.
  character(len=128) :: hostname

  integer, dimension(:), allocatable :: ele_owner
  type(halo_type) :: halo

!  call MPI_INIT ( mpi_my_error ); CHKERRQ(mpi_my_error)

! find out MY process ID, and how many processes were started.
  call MPI_COMM_RANK (MPI_COMM_WORLD, mpi_my_id, mpi_my_error); CHKERRQ(mpi_my_error)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, mpi_num_procs, mpi_my_error); CHKERRQ(mpi_my_error)

  write(mpi_num_procs_character,'(I5)') mpi_num_procs
  call hostnm(hostname)
!  write (*,*) "mpi_num_procs:",mpi_num_procs,", hostname:",trim(adjustl(hostname))

  pass = .true.

! Do everything in serial  
  if ( mpi_my_id .eq. 0 ) then  
    t_0 = mpi_wtime()
    positionsA = read_triangle_files("data/square_0_5"//"_"//trim(adjustl(mpi_num_procs_character)), dim)
    positionsB = read_triangle_files("data/square_0_9"//"_"//trim(adjustl(mpi_num_procs_character)), dim)
    serial_ele_A = ele_count(positionsA)
    serial_ele_B = ele_count(positionsB)
    mesh_serial = mpi_wtime() - t_0

!    print "(a,e25.17e3)", "Mesh input time  = ", mesh_serial
!    print "(a,i10,a,i10,a)", "Element Count (tris) A:",ele_count(positionsA),", Node Count (cells) A:",node_count(positionsA),"."
!    print "(a,i10,a,i10,a)", "Element Count (tris) B:",ele_count(positionsB),", Node Count (cells) B:",node_count(positionsB),"."

    t1 = mpi_wtime()
    do ele_A=1,ele_count(positionsA)
      tri_A%v = ele_val(positionsA, ele_A)

      do ele_B=1,ele_count(positionsB) 
        ! B. Use the libsupermesh internal triangle intersector (using derived types as input)
        tri_B%v = ele_val(positionsB, ele_B)

        call intersect_tris(tri_A, tri_B, trisC, n_trisC)

        do ele_C=1,n_trisC
          area_serial = area_serial + triangle_area(trisC(ele_C)%v)
        end do
      end do
    end do
    t2 = mpi_wtime()
    time_serial = t2 - t1

!    print "(a,e25.17e3)", "Serial intersection time  = ", time_serial
!    print "(a,e25.17e3,a)", "Serial intersection area:", area_serial,"."

    call deallocate(positionsA)
    call deallocate(positionsB)
  end if

!  ToDo TODO todo FIX HACK
!  Remove the comments when running non debug code
!  if ( mpi_my_id .eq. 0 ) then 
    allocate(areas_parallel(0:mpi_num_procs-1), times_parallel(0:mpi_num_procs-1))
!  end if
  call FLUSH()
  call MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)

!  ToDo TODO todo FIX HACK
!  Remove the DO loop and the IF check
!  do i=0,mpi_num_procs
    call FLUSH()
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
!    if ( mpi_my_id .eq. i ) then 
      area_parallel = 0.0
      write(mpi_my_id_character,'(I5)') mpi_my_id

      t_0 = mpi_wtime()
      positionsA = read_triangle_files("data/square_0_5"//"_"//trim(adjustl(mpi_num_procs_character)), dim)
      parallel_ele_A = ele_count(positionsA)


      filenameB = trim(adjustl("data/square_0_9_"))//trim(adjustl(mpi_num_procs_character))//"_"//trim(adjustl(mpi_my_id_character))
      positionsB = read_triangle_files(trim(filenameB), dim)

      call read_halo("data/square_0_9"//"_"//trim(adjustl(mpi_num_procs_character)), halo, level = 2)
      allocate(ele_owner(ele_count(positionsB)))
      call element_ownership(node_count(positionsB), reshape(positionsB%mesh%ndglno, (/ele_loc(positionsB, 1), ele_count(positionsB)/)), halo, ele_owner)
      parallel_ele_B = count(ele_owner == mpi_my_id)
      call deallocate(halo)

      mesh_parallel = mpi_wtime() - t_0

!      print "(i5,a,e25.17e3)", mpi_my_id,": Mesh input time  = ", mesh_parallel
!      print "(i5,a,i10,a,i10,a)", mpi_my_id,": Element Count (tris) A:",ele_count(positionsA),", Node Count (cells) A:",node_count(positionsA),"."
!      print "(i5,a,i10,a,i10,a)", mpi_my_id,": Element Count (tris) B:",ele_count(positionsB),", Node Count (cells) B:",node_count(positionsB),"."

      t1 = mpi_wtime()
      do ele_A=1,ele_count(positionsA)
!        write(*,*) "ele_A:",ele_A,"."
!        write(*,*) "positionsA%mesh%ndglno(ele_A) :",positionsA%mesh%ndglno(ele_A),"."
!        write(*,*) "ele_val(",ele_A,") :",ele_val(positionsA, ele_A),"."
!        write(*,*) "positionsA%val(:, ele_A) :",positionsA%val(:, :),"."
        call FLUSH()
        tri_A%v = ele_val(positionsA, ele_A)
        do ele_B=1,ele_count(positionsB) 
          ! B. Use the libsupermesh internal triangle intersector (using derived types as input)
          if(ele_owner(ele_B) /= mpi_my_id) cycle
          tri_B%v = ele_val(positionsB, ele_B)

          call intersect_tris(tri_A, tri_B, trisC, n_trisC)

          do ele_C=1,n_trisC
            area_parallel = area_parallel + triangle_area(trisC(ele_C)%v)
          end do
        end do
      end do
      t2 = mpi_wtime()
      time_parallel = t2 - t1

!      print "(i5,a,e25.17e3)", mpi_my_id, ": intersection time  = ", time_parallel
!      print "(i5,a,e25.17e3,a)", mpi_my_id, ": intersection area:", area_parallel,"."
    
      call deallocate(positionsA)
      call deallocate(positionsB)
!    end if
!    call FLUSH()
!    call MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
!  end do
  
  call FLUSH()
  call MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)

  ! Gather remote results:
  call MPI_Gather(area_parallel, 1, MPI_DOUBLE_PRECISION, &
      & areas_parallel, 1, MPI_DOUBLE_PRECISION, &
      & mpi_my_root, MPI_COMM_WORLD, mpi_my_error)
  call MPI_Gather(time_parallel, 1, MPI_DOUBLE_PRECISION, &
      & times_parallel, 1, MPI_DOUBLE_PRECISION, &
      & mpi_my_root, MPI_COMM_WORLD, mpi_my_error)
  CALL MPI_REDUCE(parallel_ele_B, local_sum, 1, MPI_INT, MPI_SUM, mpi_my_root, MPI_COMM_WORLD, mpi_my_error)

  if ( mpi_my_id .eq. 0 ) then 
!    write (*,*) times_parallel
!    write (*,*) areas_parallel
!    print "(i5,a,e25.17e3,a)", mpi_my_id, ": Total parallel intersection time:", sum(times_parallel),"."
!    print "(i5,a,e25.17e3,a)", mpi_my_id, ": Total parallel intersection area:", sum(areas_parallel),"."
    fail = fnequals(sum(areas_parallel), area_serial, tol = tol)
    call report_test("[test_parallel_partition_b areas]", fail, .false., "Should give the same areas of intersection")
    if (fail) then
      print "(a,e25.17e3,a,e25.17e3,a)", ": Total parallel intersection area:", sum(areas_parallel),&
                                       & ", total serial intersection area:",area_serial,"."
      write (*,*) areas_parallel
    end if

    fail = ( serial_ele_A .ne. parallel_ele_A )
    call report_test("[test_parallel_partition_b elementsA]", fail, .false., "Should give the same number of elements for mesh A")
    if (fail) then
      print "(a,i5,a,i5,a)", ": Total parallel elements for mesh A:",parallel_ele_A,&
                           & ", total serial elements for mesh A:",serial_ele_A,"."
    end if

    fail = ( serial_ele_B .ne. local_sum )
    call report_test("[test_parallel_partition_b elementsB]", fail, .false., "Should give the same number of elements for mesh B")
    if (fail) then
      print "(a,i5,a,i5,a)", ": Total parallel elements for mesh B:",local_sum,&
                           & ", total serial elements for mesh B:",serial_ele_B,"."
    end if
  end if
  
  deallocate(areas_parallel, times_parallel)

end subroutine test_parallel_partition_b
