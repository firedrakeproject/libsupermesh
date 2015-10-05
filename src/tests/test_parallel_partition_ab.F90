subroutine test_parallel_partition_ab

  use libsupermesh_unittest_tools
  use libsupermesh_construction
  use libsupermesh_fields
  use libsupermesh_read_triangle
  use libsupermesh_tri_intersection_module
  use libsupermesh_unittest_tools
  use libsupermesh_intersection_finder
  use libsupermesh_read_halos
  use libsupermesh_halo_ownership
  
  implicit none
  
#include <finclude/petsc.h90>

  integer :: i, j, k, nnodes, ele_A, ele_B, ele_C, n_trisC, mpi_num_procs, &
       & mpi_my_id, mpi_my_error, serial_ele_A, serial_ele_B, parallel_ele_A, &
       & parallel_ele_B, local_sum_a, local_sum_b, triangles, &
       & serial_local_iter, serial_local_iter_actual, local_iter, local_iter_actual
  integer, dimension(:), allocatable :: request
  integer, dimension(:,:), allocatable :: status
  real, dimension(:,:), allocatable :: recv_buffer
  real, dimension(:), allocatable :: send_buffer
  character(len=5) :: mpi_my_id_character
  type(tri_type) :: tri_A, tri_B
  type(tri_type), dimension(tri_buf_size) :: trisC
  logical :: pass
  real :: t_0, t1, t2, area_serial = 0.0, area_parallel = 0.0, & 
       & time_serial = 0.0, time_parallel = 0.0, mesh_serial = 0.0, mesh_parallel = 0.0
  real, dimension(:), allocatable :: areas_parallel, times_parallel
  integer, dimension(:), allocatable :: local_iters_parallel, local_iter_actual_parallel
  type(vector_field) :: positionsA, positionsB
  integer, parameter :: dim = 2
  character(len=9999) :: filenameA, filenameB
  integer, parameter :: mpi_my_root = 0
  real, parameter :: tol = 1.0e3 * epsilon(0.0)
  logical :: fail = .FALSE.
  logical, dimension(:), allocatable :: partition_intersection
  real, dimension(:,:), allocatable :: bbox_a, bbox_b
  real, dimension(:,:,:), allocatable :: parallel_bbox_a, parallel_bbox_b
  integer, dimension(:), allocatable  :: parallel_ele_B_array

  integer, dimension(:), allocatable :: ele_ownerA, ele_ownerB
  type(halo_type) :: halo

  type(intersections), dimension(:), allocatable :: map_AB

!  call MPI_INIT ( mpi_my_error ); CHKERRQ(mpi_my_error)

! find out MY process ID, and how many processes were started.
  CALL MPI_COMM_RANK (MPI_COMM_WORLD, mpi_my_id, mpi_my_error); CHKERRQ(mpi_my_error)
  CALL MPI_COMM_SIZE (MPI_COMM_WORLD, mpi_num_procs, mpi_my_error); CHKERRQ(mpi_my_error)

  allocate(request(2*mpi_num_procs))
  allocate(status(MPI_STATUS_SIZE,2*mpi_num_procs))
  allocate(partition_intersection(0:mpi_num_procs-1))
  partition_intersection = .FALSE.
  pass = .TRUE.
  serial_local_iter = 0
  serial_local_iter_actual = 0
  local_iter = 0
  local_iter_actual = 0
  k = 0
  triangles = 0

! Do everything in serial  
  if ( mpi_my_id .eq. 0 ) then
    t_0 = mpi_wtime()
    positionsA = read_triangle_files("data/square_0_2", dim)
!    positionsB = read_triangle_files("data/square_0_1", dim)
    positionsB = read_triangle_files("data/square_0_002", dim)
    serial_ele_A = ele_count(positionsA)
    serial_ele_B = ele_count(positionsB)
    mesh_serial = mpi_wtime() - t_0

!    print "(a,e25.17e3)", "Mesh input time  = ", mesh_serial
!    print "(a,i10,a,i10,a)", "Element Count (tris) A:",ele_count(positionsA),", Node Count (cells) A:",node_count(positionsA),"."
!    print "(a,i10,a,i10,a)", "Element Count (tris) B:",ele_count(positionsB),", Node Count (cells) B:",node_count(positionsB),"."
    
    t1 = mpi_wtime()
    allocate(map_AB(serial_ele_A))
    call intersection_finder(positionsA%val, reshape(positionsA%mesh%ndglno, (/ele_loc(positionsA, 1), serial_ele_A/)), &
                           & positionsB%val, reshape(positionsB%mesh%ndglno, (/ele_loc(positionsB, 1), serial_ele_B/)), &
                           & map_AB)

    do ele_A=1,ele_count(positionsA)
      tri_A%v = ele_val(positionsA, ele_A)

      do i = 1, map_AB(ele_A)%n
        ele_B = map_AB(ele_A)%v(i)
        ! B. Use the libSuperMesh internal triangle intersector (using derived types as input)
        tri_B%v = ele_val(positionsB, ele_B)

        call intersect_tris(tri_A, tri_B, trisC, n_trisC)

        serial_local_iter = serial_local_iter + 1

        do ele_C=1,n_trisC
          area_serial = area_serial + triangle_area(trisC(ele_C)%v)
          serial_local_iter_actual = serial_local_iter_actual + 1
        end do
      end do
    end do
    
    deallocate(map_AB)
    t2 = mpi_wtime()
    time_serial = t2 - t1

!    print "(a,e25.17e3)", "Serial intersection time  = ", time_serial
!    print "(a,e25.17e3,a)", "Serial intersection area:", area_serial,"."

    call deallocate(positionsA)
    call deallocate(positionsB)
  end if

  allocate(parallel_ele_B_array(0:mpi_num_procs-1))
!  ToDo TODO todo FIX HACK
!  Remove the comments when running non debug code
!  if ( mpi_my_id .eq. 0 ) then 
    allocate(areas_parallel(0:mpi_num_procs-1), times_parallel(0:mpi_num_procs-1))
    allocate(local_iters_parallel(0:mpi_num_procs-1))
    allocate(local_iter_actual_parallel(0:mpi_num_procs-1))
!  end if
  CALL FLUSH()
  CALL MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
  
!  ToDo TODO todo FIX HACK
!  Remove the DO loop and the IF check
!  do i=0,mpi_num_procs
!    call FLUSH()
!    call MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
!    if ( mpi_my_id .eq. i ) then
      area_parallel = 0.0
      write(mpi_my_id_character,'(I5)') mpi_my_id

      t_0 = mpi_wtime()
      filenameA = trim(adjustl("data/square_0_2_"))//trim(adjustl(mpi_my_id_character))

!      print "(a,i5,a,i10,a,a,a)", "MPI Process:",mpi_my_id,", has nodes:",nnodes,", read from:",trim(filenameA),"."
      CALL FLUSH()

      positionsA = read_triangle_files(trim(filenameA), dim)
      parallel_ele_A = ele_count(positionsA)

      call read_halo("data/square_0_2", halo, level = 2)
      allocate(ele_ownerA(ele_count(positionsA)))
      call element_ownership(node_count(positionsA), reshape(positionsA%mesh%ndglno, (/ele_loc(positionsA, 1), ele_count(positionsA)/)), halo, ele_ownerA)
      parallel_ele_A = count(ele_ownerA == mpi_my_id)
      call deallocate(halo)


!      filenameB = trim(adjustl("data/square_0_1_"))//trim(adjustl(mpi_my_id_character))
      filenameB = trim(adjustl("data/square_0_002_"))//trim(adjustl(mpi_my_id_character))
      positionsB = read_triangle_files(trim(filenameB), dim)

!      call read_halo("data/square_0_1", halo, level = 2)
      call read_halo("data/square_0_002", halo, level = 2)
      allocate(ele_ownerB(ele_count(positionsB)))
      call element_ownership(node_count(positionsB), reshape(positionsB%mesh%ndglno, (/ele_loc(positionsB, 1), ele_count(positionsB)/)), halo, ele_ownerB)
      parallel_ele_B = count(ele_ownerB == mpi_my_id)
      call deallocate(halo)

      mesh_parallel = mpi_wtime() - t_0

!      print "(i5,a,e25.17e3)", mpi_my_id,": Mesh input time  = ", mesh_parallel
!      print "(i5,a,i10,a,i10,a)", mpi_my_id,": Element Count (tris) A:",ele_count(positionsA),", Node Count (cells) A:",node_count(positionsA),"."
!      print "(i5,a,i10,a,i10,a)", mpi_my_id,": Element Count (tris) B:",ele_count(positionsB),", Node Count (cells) B:",node_count(positionsB),"."

      t1 = mpi_wtime()
      allocate(bbox_a(positionsA%dim,2))
      bbox_a = partition_bbox(positionsA, ele_ownerA, mpi_my_id)
      allocate(bbox_b(positionsB%dim,2))
      bbox_b = partition_bbox(positionsB, ele_ownerB, mpi_my_id)
!      write (*,*) "                       Xmin,                 Xmax,                   Ymin,             Ymax"
!      print "(i5,a,F19.17,2x,F19.17,4x,F19.17,2x,F19.17,a)", mpi_my_id," bbox_a    :",bbox_a,"."
!      print "(i5,a,F19.17,2x,F19.17,4x,F19.17,2x,F19.17,a)", mpi_my_id," bbox_b    :",bbox_b,"."
!    end if
!    CALL FLUSH()
!    CALL MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
!  end do

  allocate(parallel_bbox_a(2, positionsA%dim, 0:mpi_num_procs-1))
  allocate(parallel_bbox_b(2, positionsB%dim, 0:mpi_num_procs-1))

  CALL MPI_Allgather(bbox_a, 4, MPI_DOUBLE_PRECISION, &
    &    parallel_bbox_a, 4, MPI_DOUBLE_PRECISION,    &
    &    MPI_COMM_WORLD, mpi_my_error)

  CALL MPI_Allgather(bbox_b, 4, MPI_DOUBLE_PRECISION, &
    &    parallel_bbox_b, 4, MPI_DOUBLE_PRECISION,    &
    &    MPI_COMM_WORLD, mpi_my_error)

  CALL MPI_Allgather(parallel_ele_B, 1, MPI_INTEGER, &
    &    parallel_ele_B_array, 1, MPI_INTEGER,    &
    &    MPI_COMM_WORLD, mpi_my_error)
    
  ! JM: Work out which elements to send using the bounding boxes

                       ! JM: Per-processing sizing required
  allocate(recv_buffer(maxval(parallel_ele_B_array)*(positionsB%dim + 1)*2,0:mpi_num_procs-1))
  allocate(send_buffer(parallel_ele_B*(positionsB%dim + 1)*2))
  recv_buffer = -9.0
  send_buffer = -100.0
  i=1
  do ele_B=1,ele_count(positionsB) 
    if(ele_ownerB(ele_B) /= mpi_my_id) cycle
    ! JM: Somewhere in the send loop
    !if(.not. bboxes_intersect(ele_val(positionsB, ele_B), parallel_bbox_a(:, :, proc))) cycle
    
    tri_B%v = ele_val(positionsB, ele_B)
    send_buffer(i) = tri_B%v(1,1)
    i = i + 1
    send_buffer(i) = tri_B%v(2,1)
    i = i + 1
    send_buffer(i) = tri_B%v(1,2)
    i = i + 1
    send_buffer(i) = tri_B%v(2,2)
    i = i + 1
    send_buffer(i) = tri_B%v(1,3)
    i = i + 1
    send_buffer(i) = tri_B%v(2,3)
    i = i + 1
  end do

!  ToDo TODO todo FIX HACK
!  Remove the DO loop and the IF check
!  do i=0,mpi_num_procs-1
!    CALL FLUSH()
!    CALL MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
!    if ( mpi_my_id .eq. i ) then
!      write (*,*) ""
!      write(*,*) mpi_my_id," parallel_bbox_a   :",parallel_bbox_a,"."
!      print "(i5,a,F19.17,2x,F19.17,8x,F19.17,2x,F19.17,8x,F19.17,2x,F19.17,8x,F19.17,2x,F19.17,8x,F19.17,2x,F19.17,8x,F19.17,2x,F19.17,8x,F19.17,2x,F19.17,8x,F19.17,2x,F19.17,a)", mpi_my_id," parallel_bbox_a :",parallel_bbox_a,"."
!      write(*,*) mpi_my_id," parallel_bbox_b   :",parallel_bbox_b,"."
!      print "(i5,a,F19.17,2x,F19.17,8x,F19.17,2x,F19.17,8x,F19.17,2x,F19.17,8x,F19.17,2x,F19.17,8x,F19.17,2x,F19.17,8x,F19.17,2x,F19.17,8x,F19.17,2x,F19.17,8x,F19.17,2x,F19.17,a)", mpi_my_id," parallel_bbox_b :",parallel_bbox_b,"."

      ! For the time being copy Partition B to all Processes
      do j=0,mpi_num_procs-1
        if (mpi_my_id .eq. j) then
          recv_buffer(:,j) = 2.0
          cycle
        end if
      end do

      ! Check bounding box intersection
      do j=0,mpi_num_procs-1
        if (mpi_my_id .eq. j) cycle
        if (bboxes_intersect(bbox_a, parallel_bbox_b(:, :, j))) then
!         write (*,*) mpi_my_id," My partition of meshA intersects with ",j," partition of meshB."
!          print "(i5,a,i2,a)", mpi_my_id,", RECVing from ",j,"."
          k = k + 1
! ! !          print "(i5,a,i5,a,i5,a)", mpi_my_id," j:",j,",k:",k,"."
          CALL MPI_Irecv(recv_buffer(1,j), & 
               & parallel_ele_B_array(j) * (positionsB%dim + 1)*2, &
               & MPI_DOUBLE_PRECISION, &
               & j, MPI_ANY_TAG, MPI_COMM_WORLD, &
               & request(k), mpi_my_error)
          partition_intersection(j) = .TRUE.
        end if

        if (bboxes_intersect(parallel_bbox_a(:, :, j), bbox_b)) then
!          print "(i5,a,i2,a)", mpi_my_id,", SENDing to ",j,"."
          k = k + 1
! ! !          print "(i5,a,i5,a,i5,a)", mpi_my_id," j:",j,",k:",k,"."
          CALL MPI_Isend(send_buffer, &
! ! !             & parallel_ele_B_array(i) * (positionsB%dim + 1), &
             & parallel_ele_B * (positionsB%dim + 1)*2, &
             & MPI_DOUBLE_PRECISION, &
             & j, 0, MPI_COMM_WORLD, request(k), mpi_my_error)
        end if

! ! !          print "(i5,a,i5,a,i5,a,i5,a)", mpi_my_id," j:",j,",k:",k,", parallel_ele_B_array(j):",parallel_ele_B_array(j),"."
!           k = k + 1
! ! !          print "(i5,a,i5,a,i5,a)", mpi_my_id," j:",j,",k:",k,"."
!           CALL MPI_Irecv(recv_buffer(1,j), & 
!              & parallel_ele_B_array(j) * (positionsB%dim + 1)*2, &
!              & MPI_DOUBLE_PRECISION, &
!              & j, MPI_ANY_TAG, MPI_COMM_WORLD, &
!              & request(k), mpi_my_error)
!           k = k + 1
! ! !          print "(i5,a,i5,a,i5,a)", mpi_my_id," j:",j,",k:",k,"."
!           CALL MPI_Isend(send_buffer, &
! ! !             & parallel_ele_B_array(i) * (positionsB%dim + 1), &
!              & parallel_ele_B * (positionsB%dim + 1)*2, &
!              & MPI_DOUBLE_PRECISION, &
!              & j, 0, MPI_COMM_WORLD, request(k), mpi_my_error)
      end do
!      print "(i5,a,F19.17,a,F19.17,a,F19.17,a,F19.17,a)", mpi_my_id,", bbox_a(1,1):",bbox_a(1,1),", bbox_a(2,1):",bbox_a(2,1),", bbox_a(1,2):",bbox_a(1,2),", bbox_a(2,2):",bbox_a(2,2),"."
!    end if
!    CALL FLUSH()
!    call MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
!  end do

!  CALL FLUSH()
!  CALL MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
!  print "(i5,a,i5,a)", mpi_my_id,", k:",k,"."
!  CALL MPI_Wait(send_request(1:k), status(1:k), mpi_my_error);
!  CALL MPI_Wait(recv_request(1:k), status(1:k), mpi_my_error);
  CALL MPI_Waitall(k, request(1:k), status(:,1:k), mpi_my_error)
!  CALL MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)

!  ToDo TODO todo FIX HACK
!  Remove the DO loop and the IF check
!  do i=0,mpi_num_procs
!    CALL FLUSH()
!    CALL MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
!    if ( mpi_my_id .eq. i ) then
!      print "(i5,a,i10,a,i10,a,i10,a)", mpi_my_id,": Total Element Count (tris) A:",ele_count(positionsA),", My Element Count (tris) A:",count(ele_ownerA == mpi_my_id),", Total Node Count (cells) A:",node_count(positionsA),"."
!      print "(i5,a,i10,a,i10,a,i10,a)", mpi_my_id,": Total Element Count (tris) B:",ele_count(positionsB),", My Element Count (tris) B:",count(ele_ownerB == mpi_my_id),", Total Node Count (cells) B:",node_count(positionsB),"."
!      write(*,*) mpi_my_id," , ",parallel_ele_B," , ",parallel_ele_B_array
!      CALL FLUSH()
!if (i == 0 .OR. i == 1) then
!if(i == 0) then
!      do j=0,mpi_num_procs-1
!        write(*,*) mpi_my_id," ", j , "."
!        if (j == 0) then 
!          write(*,*) mpi_my_id,"  send_buffer:"
!          write(*,*) send_buffer
!        end if
!        write(*,*) mpi_my_id,"  recv_buffer:"
!        write(*,*) recv_buffer(:,j)
!        write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!      end do
!!!      call EXIT(-1)
!end if
      do ele_A=1,ele_count(positionsA)
        if(ele_ownerA(ele_A) /= mpi_my_id) then
!          print *, mpi_my_id,": Element:",ele_A,", does NOT belong to me. Element values:",ele_val(positionsA, ele_A),"."
          cycle
        end if
!        if(ele_ownerA(ele_A) == mpi_my_id) then
!          print *, mpi_my_id,": Element:",ele_A,", BELONGS to me. Element values:",ele_val(positionsA, ele_A),"."
!        end if
!        write(*,*) "ele_A:",ele_A,"."
!        write(*,*) "positionsA%mesh%ndglno(ele_A) :",positionsA%mesh%ndglno(ele_A),"."
!        print "(i2,a,i3,a,F19.17,2x,F19.17,8x,F19.17,2x,F19.17,8x,F19.17,2x,F19.17,a)", mpi_my_id," ele_val(",ele_A,"):",ele_val(positionsA, ele_A),"."
!        write(*,*) "ele_val(",ele_A,") :",ele_val(positionsA, ele_A),"."
!        write(*,*) "positionsA%val(:, ele_A) :",positionsA%val(:, :),"."
        CALL FLUSH()
        tri_A%v = ele_val(positionsA, ele_A)
        do ele_B=1,ele_count(positionsB) 
          if(ele_ownerB(ele_B) /= mpi_my_id) cycle
          ! B. Use the libSuperMesh internal triangle intersector (using derived types as input)
          tri_B%v = ele_val(positionsB, ele_B)

          local_iter = local_iter + 1

          call intersect_tris(tri_A, tri_B, trisC, n_trisC)

          do ele_C=1,n_trisC
            area_parallel = area_parallel + triangle_area(trisC(ele_C)%v)
          end do
        end do

!        print "(i5,a,i5,a,F19.17,a)", mpi_my_id,", i:",i,", area_parallel (TEMP):",area_parallel,"."

        do j=0,mpi_num_procs-1
          if ( j .eq. i ) cycle
          if ( partition_intersection(j) .eqv. .FALSE. ) cycle
!          print "(i5,a,i5,a)", mpi_my_id,", Using buffer RECVed from :",j,"."
          do k=1,parallel_ele_B_array(j) * (positionsB%dim + 1)*2,6
             tri_B%v(1,1) = recv_buffer(k,j)
             tri_B%v(2,1) = recv_buffer(k+1,j)
             tri_B%v(1,2) = recv_buffer(k+2,j)
             tri_B%v(2,2) = recv_buffer(k+3,j)
             tri_B%v(1,3) = recv_buffer(k+4,j)
             tri_B%v(2,3) = recv_buffer(k+5,j)

             local_iter = local_iter + 1
!             write (*,*) mpi_my_id,", tri_B%v:",tri_B%v

             call intersect_tris(tri_A, tri_B, trisC, n_trisC)

             do ele_C=1,n_trisC
              area_parallel = area_parallel + triangle_area(trisC(ele_C)%v)
              local_iter_actual = local_iter_actual + 1
            end do
!            print "(i5,a,i5,a,i5,a,i5,a,F19.17,a)", mpi_my_id,", i:",i,", j:",j,", k:",k,", area_parallel:",area_parallel,"."
          end do
!          print "(i5,a,i5,a,i5,a)", mpi_my_id,", parallel_ele_B_array(j):",parallel_ele_B_array(j),"."
!          print "(i5,a,i5,a,i5,a)", mpi_my_id,", parallel_ele_B_array(j) * (positionsB%dim + 1)*2:",parallel_ele_B_array(j) * (positionsB%dim + 1)*2,"."
          triangles = triangles + k/6
!          print "(i5,a,i5,a,i5,a)", mpi_my_id," ele_A:",ele_A,", done ",triangles,"."
!          call EXIT(-1)

        end do
      end do
      t2 = mpi_wtime()
      time_parallel = t2 - t1

!      print "(i5,a,e25.17e3)", mpi_my_id, ": intersection time  = ", time_parallel
!      print "(i5,a,e25.17e3,a)", mpi_my_id, ": intersection area:", area_parallel,"."

      call deallocate(positionsA)
      call deallocate(positionsB)
      deallocate(ele_ownerA)
      deallocate(ele_ownerB)
!    end if
!    call FLUSH()
!    call MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
!  end do

!  print "(i5,a,i5,a)", mpi_my_id," done ",triangles,"."
!  CALL FLUSH()
!  CALL MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
!call exit(-1)
  ! Gather remote results:
  CALL MPI_Gather(area_parallel, 1, MPI_DOUBLE_PRECISION, &
      & areas_parallel, 1, MPI_DOUBLE_PRECISION, &
      & mpi_my_root, MPI_COMM_WORLD, mpi_my_error)

  CALL MPI_Gather(time_parallel, 1, MPI_DOUBLE_PRECISION, &
      & times_parallel, 1, MPI_DOUBLE_PRECISION, &
      & mpi_my_root, MPI_COMM_WORLD, mpi_my_error)

  CALL MPI_Gather(local_iter, 1, MPI_INTEGER, &
      & local_iters_parallel, 1, MPI_INTEGER, &
      & mpi_my_root, MPI_COMM_WORLD, mpi_my_error)

  CALL MPI_Gather(local_iter_actual, 1, MPI_INTEGER, &
      & local_iter_actual_parallel, 1, MPI_INTEGER, &
      & mpi_my_root, MPI_COMM_WORLD, mpi_my_error)

  CALL MPI_REDUCE(parallel_ele_A, local_sum_a, 1, MPI_INT, MPI_SUM, mpi_my_root, MPI_COMM_WORLD, mpi_my_error)

  CALL MPI_REDUCE(parallel_ele_B, local_sum_b, 1, MPI_INT, MPI_SUM, mpi_my_root, MPI_COMM_WORLD, mpi_my_error)

  if ( mpi_my_id .eq. 0 ) then 
!    write (*,*) times_parallel
!    write (*,*) areas_parallel
    print "(i5,a,i15,a,i15,a)", mpi_my_id, ": Ele A:",serial_ele_A,", Ele B:",serial_ele_B,"."
    print "(i5,a,e25.17e3,a)", mpi_my_id, ": Total serial intersection time  :", time_serial,"."
    print "(i5,a,e25.17e3,a)", mpi_my_id, ": Total parallel intersection time:", sum(times_parallel),"."
    print "(i5,a,e25.17e3,a)", mpi_my_id, ": Max parallel intersection time  :", maxval(times_parallel),"."
    print "(i5,a,e25.17e3,a)", mpi_my_id, ": Total serial intersection area  :", area_serial,"."
    print "(i5,a,e25.17e3,a)", mpi_my_id, ": Total parallel intersection area:", sum(areas_parallel),"."
!    print "(i5,a,e25.17e3,a)", mpi_my_id, ": Total parallel intersection area:", sum(areas_parallel),"."
    print "(i5,a,i15,a,i15,a)", mpi_my_id, ": Serial iters         :",serial_local_iter,"."
    print "(i5,a,i15,a,i15,a)", mpi_my_id, ": Parallel iters       :",sum(local_iters_parallel),"."
    print "(i5,a,i15,a,i15,a)", mpi_my_id, ": Serial actual iters  :",serial_local_iter_actual,"."
    print "(i5,a,i15,a,i15,a)", mpi_my_id, ": Parallel actual iters:",sum(local_iter_actual_parallel),"."
    fail = fnequals(sum(areas_parallel), area_serial, tol = tol)
    call report_test("[test_parallel_partition_ab areas]", fail, .FALSE., "Should give the same areas of intersection")
    if (fail) then
      print "(a,e25.17e3,a,e25.17e3,a)", ": Total parallel intersection area:", sum(areas_parallel),&
                                       & ", total serial intersection area:",area_serial,"."
      write (*,*) "areas_parallel:",areas_parallel
      write (*,*) "local_iters_parallel:",local_iters_parallel
      print "(a,i15,a,i15,a)", "Serial iters:",serial_local_iter,", parallel iters:",sum(local_iters_parallel),"."
    end if

    fail = ( serial_ele_A .ne. local_sum_a )
    call report_test("[test_parallel_partition_ab elementsA]", fail, .FALSE., "Should give the same number of elements for mesh A")
    if (fail) then
      print "(a,i5,a,i5,a)", ": Total parallel elements for mesh A:",local_sum_a,&
                           & ", total serial elements for mesh A:",serial_ele_A,"."
    end if

    fail = ( serial_ele_B .ne. local_sum_b )
    call report_test("[test_parallel_partition_ab elementsB]", fail, .FALSE., "Should give the same number of elements for mesh B")
    if (fail) then
      print "(a,i10,a,i10,a)", ": Total parallel elements for mesh B:",local_sum_b,&
                           & ", total serial elements for mesh B:",serial_ele_B,"."
    end if
  end if

  deallocate(areas_parallel, times_parallel)
  deallocate(parallel_bbox_a, parallel_bbox_b)
  deallocate(recv_buffer, send_buffer, parallel_ele_B_array)
  deallocate(request, status, partition_intersection)

  call cintersection_finder_reset(nnodes)

end subroutine test_parallel_partition_ab
