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

  type pointer_real
    real, dimension(:), pointer :: p
  end type pointer_real

  type pointer_integer
    integer, dimension(:), pointer :: p
  end type pointer_integer

#include <finclude/petsc.h90>

  integer :: i, j, k, l, nnodes, ele_A, ele_B, ele_C, n_trisC, mpi_num_procs, &
       & mpi_my_id, mpi_my_error, serial_ele_A, serial_ele_B, parallel_ele_A, &
       & parallel_ele_B
  integer :: local_sum_a, local_sum_b, triangles, &
       & serial_local_iter, serial_local_iter_actual, local_iter, local_iter_actual
  logical :: pass
  real :: t0, t1, t2, area_serial = 0.0, area_parallel = 0.0, & 
       & time_serial = 0.0, local_time = 0.0, mesh_serial = 0.0, &
       & local_time_intersection_only = 0.0, local_mesh = 0.0, local_other_time = 0.0
  real, dimension(:), allocatable :: areas_parallel, times_parallel, times_intersection_only_parallel, other_times_parallel
  integer, dimension(:), allocatable :: iters_parallel, iter_actual_parallel
  real, dimension(:), allocatable :: mesh_parallel

  integer, dimension(:), allocatable   :: number_of_elements_to_receive
  integer, dimension(:), allocatable   :: temp_elements_uns
  type(pointer_integer), dimension(:), allocatable  :: elements_uns
!  integer, dimension(:,:), allocatable :: elements_uns
  integer, dimension(:), allocatable   :: request
  integer, dimension(:,:), allocatable :: status
!  real, dimension(:,:), allocatable    :: recv_buffer
  type(pointer_real), dimension(:), allocatable  :: recv_buffer, send_buffer
  real, dimension(:), allocatable  :: temp_buffer
  character(len=5) :: mpi_my_id_character, mpi_num_procs_character
  type(tri_type) :: tri_A, tri_B
  type(tri_type), dimension(tri_buf_size) :: trisC
  type(vector_field) :: positionsA, positionsB
  integer, parameter :: dim = 2
  character(len=9999) :: filenameA, filenameB
  integer, parameter :: mpi_my_root = 0
  real, parameter :: tol = 1.0e3 * epsilon(0.0)
  logical :: fail = .FALSE.
  logical :: toprint = .TRUE.
  logical, dimension(:), allocatable :: partition_intersection
  real, dimension(:,:), allocatable :: bbox_a, bbox_b
  real, dimension(:,:,:), allocatable :: parallel_bbox_a, parallel_bbox_b
  integer, dimension(:), allocatable  :: parallel_ele_B_array

  integer, dimension(:), allocatable :: ele_ownerA, ele_ownerB
  type(halo_type) :: halo

  type(intersections), dimension(:), allocatable :: map_AB
  character(len=128) :: hostname

!  call MPI_INIT ( mpi_my_error ); CHKERRQ(mpi_my_error)

! find out MY process ID, and how many processes were started.
  CALL MPI_COMM_RANK (MPI_COMM_WORLD, mpi_my_id, mpi_my_error); CHKERRQ(mpi_my_error)
  CALL MPI_COMM_SIZE (MPI_COMM_WORLD, mpi_num_procs, mpi_my_error); CHKERRQ(mpi_my_error)

  write(mpi_num_procs_character,'(I5)') mpi_num_procs
  call hostnm(hostname)
!  write (*,*) "mpi_num_procs:",mpi_num_procs,", hostname:",trim(adjustl(hostname))

  allocate(number_of_elements_to_receive(0:mpi_num_procs-1))
  number_of_elements_to_receive = 0
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
    t0 = mpi_wtime()
    positionsA = read_triangle_files("data/square_0_2"//"_"//trim(adjustl(mpi_num_procs_character)), dim)
!    positionsA = read_triangle_files("data/square_0_02"//"_"//trim(adjustl(mpi_num_procs_character)), dim)
    positionsB = read_triangle_files("data/square_0_1"//"_"//trim(adjustl(mpi_num_procs_character)), dim)
!    positionsB = read_triangle_files("data/square_0_004"//"_"//trim(adjustl(mpi_num_procs_character)), dim)
!    positionsB = read_triangle_files("data/square_0_002"//"_"//trim(adjustl(mpi_num_procs_character)), dim)
    serial_ele_A = ele_count(positionsA)
    serial_ele_B = ele_count(positionsB)
    mesh_serial = mpi_wtime() - t0

!    print "(a,e25.17e3)", "Mesh input time  = ", mesh_serial
!    print "(a,i10,a,i10,a)", "Element Count (tris) A:",ele_count(positionsA),", Node Count (cells) A:",node_count(positionsA),"."
!    print "(a,i10,a,i10,a)", "Element Count (tris) B:",ele_count(positionsB),", Node Count (cells) B:",node_count(positionsB),"."

    t1 = mpi_wtime()
!    Use the intersection finder!!
!    allocate(map_AB(serial_ele_A))
!    call intersection_finder(positionsA%val, reshape(positionsA%mesh%ndglno, (/ele_loc(positionsA, 1), serial_ele_A/)), &
!                           & positionsB%val, reshape(positionsB%mesh%ndglno, (/ele_loc(positionsB, 1), serial_ele_B/)), &
!                           & map_AB)

    do ele_A=1,ele_count(positionsA)
      tri_A%v = ele_val(positionsA, ele_A)
!      if ( ele_A < 5 ) write(*,*) "A: ele_A:",ele_A,", ",tri_A%v

      do ele_B=1,ele_count(positionsB)
!    Use the intersection finder!!
!      do i = 1, map_AB(ele_A)%n
!        ele_B = map_AB(ele_A)%v(i)
        ! B. Use the libSuperMesh internal triangle intersector (using derived types as input)
        tri_B%v = ele_val(positionsB, ele_B)

        call intersect_tris(tri_A, tri_B, trisC, n_trisC)

        serial_local_iter = serial_local_iter + 1

        do ele_C=1,n_trisC
          area_serial = area_serial + triangle_area(trisC(ele_C)%v)
          serial_local_iter_actual = serial_local_iter_actual + 1
!          if ( ele_A < 5 ) then
!            write(*,*) "B: ele_B:",ele_B,", ",tri_B%v,", n_trisC:",n_trisC
!          end if
        end do
      end do
!      if ( ele_A < 7 ) write (*,*) ""
!      print "(a,i5,a,F19.17,a)", "   ele_A:",ele_A,", area_serial:",area_serial,"."
    end do

!    Use the intersection finder!!
!    deallocate(map_AB)
    t2 = mpi_wtime()
    time_serial = t2 - t1

!    print "(a,e25.17e3)", "Serial intersection time  = ", time_serial
!    print "(a,e25.17e3,a)", "Serial intersection area:", area_serial,"."

    call deallocate(positionsA)
    call deallocate(positionsB)
  end if

  CALL FLUSH()
  CALL MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)

  allocate(parallel_ele_B_array(0:mpi_num_procs-1))
!  ToDo TODO todo FIX HACK
!  Remove the comments when running non debug code
!  if ( mpi_my_id .eq. 0 ) then 
    allocate(areas_parallel(0:mpi_num_procs-1), times_parallel(0:mpi_num_procs-1), other_times_parallel(0:mpi_num_procs-1))
    allocate(iters_parallel(0:mpi_num_procs-1))
    allocate(iter_actual_parallel(0:mpi_num_procs-1))
    allocate(mesh_parallel(0:mpi_num_procs-1))
    allocate(times_intersection_only_parallel(0:mpi_num_procs-1))
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

      t0 = mpi_wtime()
      filenameA = trim(adjustl("data/square_0_2_"))//trim(adjustl(mpi_num_procs_character))//"_"//trim(adjustl(mpi_my_id_character))
!      filenameA = trim(adjustl("data/square_0_02_"))//trim(adjustl(mpi_num_procs_character))//"_"//trim(adjustl(mpi_my_id_character))

!      print "(a,i5,a,i10,a,a,a)", "MPI Process:",mpi_my_id,", has nodes:",nnodes,", read from:",trim(filenameA),"."
      CALL FLUSH()

      positionsA = read_triangle_files(trim(filenameA), dim)
      parallel_ele_A = ele_count(positionsA)

      call read_halo("data/square_0_2"//"_"//trim(adjustl(mpi_num_procs_character)), halo, level = 2)
!      call read_halo("data/square_0_02"//"_"//trim(adjustl(mpi_num_procs_character)), halo, level = 2)
      allocate(ele_ownerA(ele_count(positionsA)))
      call element_ownership(node_count(positionsA), reshape(positionsA%mesh%ndglno, (/ele_loc(positionsA, 1), ele_count(positionsA)/)), halo, ele_ownerA)
      parallel_ele_A = count(ele_ownerA == mpi_my_id)
      call deallocate(halo)


      filenameB = trim(adjustl("data/square_0_1_"))//trim(adjustl(mpi_num_procs_character))//"_"//trim(adjustl(mpi_my_id_character))
!      filenameB = trim(adjustl("data/square_0_004_"))//trim(adjustl(mpi_num_procs_character))//"_"//trim(adjustl(mpi_my_id_character))
      positionsB = read_triangle_files(trim(filenameB), dim)

      call read_halo("data/square_0_1"//"_"//trim(adjustl(mpi_num_procs_character)), halo, level = 2)
!      call read_halo("data/square_0_004"//"_"//trim(adjustl(mpi_num_procs_character)), halo, level = 2)
      allocate(ele_ownerB(ele_count(positionsB)))
      call element_ownership(node_count(positionsB), reshape(positionsB%mesh%ndglno, (/ele_loc(positionsB, 1), ele_count(positionsB)/)), halo, ele_ownerB)
      parallel_ele_B = count(ele_ownerB == mpi_my_id)
      call deallocate(halo)

      local_mesh = mpi_wtime() - t0

!      print "(i5,a,e25.17e3)", mpi_my_id,": Mesh input time  = ", mesh_parallel
!      print "(i5,a,i10,a,i10,a)", mpi_my_id,": Element Count (tris) A:",ele_count(positionsA),", Node Count (cells) A:",node_count(positionsA),"."
!      print "(i5,a,i10,a,i10,a)", mpi_my_id,": Element Count (tris) B:",ele_count(positionsB),", Node Count (cells) B:",node_count(positionsB),"."

      t1 = mpi_wtime()
      allocate(bbox_a(positionsA%dim,2))
      bbox_a = partition_bbox(positionsA, ele_ownerA, mpi_my_id)
      allocate(bbox_b(positionsB%dim,2))
      bbox_b = partition_bbox(positionsB, ele_ownerB, mpi_my_id)
!      write (*,*) mpi_my_id," bbox_a:",bbox_a
!      write (*,*) mpi_my_id," bbox_b:",bbox_b
!    end if
!    CALL FLUSH()
!    CALL MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
!  end do

  allocate(parallel_bbox_a(2, positionsA%dim, 0:mpi_num_procs-1))
  allocate(parallel_bbox_b(2, positionsB%dim, 0:mpi_num_procs-1))

  ! Required GLOBAL communication
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
!  allocate(elements_uns(maxval(parallel_ele_B_array)*(positionsB%dim + 1)*2,0:mpi_num_procs-1))
  allocate(elements_uns(0:mpi_num_procs-1))
  do i=0,size(elements_uns(:))-1
    nullify(elements_uns(i)%p)
  end do

  allocate(send_buffer(0:mpi_num_procs-1))
  do i=0,size(send_buffer(:))-1
    nullify(send_buffer(i)%p)
  end do

  allocate(recv_buffer(0:mpi_num_procs-1))
  do i=0,size(recv_buffer(:))-1
    nullify(recv_buffer(i)%p)
  end do

  k = 0
!  ToDo TODO todo FIX HACK
!  Remove the DO loop and the IF check
!  do i=0,mpi_num_procs-1
!    CALL FLUSH()
!    CALL MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
!    if ( mpi_my_id .eq. i ) then
      do j=0,mpi_num_procs-1
        l = 0

        if (bboxes_intersect(bbox_a, parallel_bbox_b(:, :, j))) then
          if ( j /= mpi_my_id ) then
            k = k + 1
            partition_intersection(j) = .TRUE.
! ! ! !            print "(i5,a,i2,a)", mpi_my_id,", we are RECVing from ",j,"."

            CALL MPI_Irecv(number_of_elements_to_receive(j), 1, &
               & MPI_INTEGER, &
               & j, MPI_ANY_TAG, MPI_COMM_WORLD, &
               & request(k), mpi_my_error)
          end if
        end if

!        if ((i .eq. 2) .and. (j .eq. 3) ) write (*,*) mpi_my_id,", j:",j,", parallel_bbox_a(:, :, j):",parallel_bbox_a(:, :, j)
        if (bboxes_intersect(parallel_bbox_a(:, :, j), bbox_b)) then
          if ( j /= mpi_my_id ) then
            k = k + 1
!            print "(i5,a,i2,a,i15,a)", mpi_my_id,", we are SENDing to ",j,". We have ",parallel_ele_B," elements in total."
            allocate(temp_elements_uns(ele_count(positionsB)))
            do ele_B=1,ele_count(positionsB) 
              if(ele_ownerB(ele_B) /= mpi_my_id) cycle
!              if ((i .eq. 2) .and. (j .eq. 3) )  print "(i5,a,i5,a,i5,a,F19.17,2x,F19.17,8x,F19.17,2x,F19.17,8x,F19.17,2x,F19.17)", mpi_my_id," TESTING (ele_B:",ele_B," for process:",j,", element:",ele_val(positionsB, ele_B)
!              if ( (i .eq. 2 ) .and. (j .eq. 3) .and. ((ele_B .eq. 10 ) .OR. (ele_B .eq. 9 )) ) write (*,*) mpi_my_id," (ele_B:",ele_B," bbox:",bbox(ele_val(positionsB, ele_B))
              if(.not. bboxes_intersect(bbox(ele_val(positionsB, ele_B)), parallel_bbox_a(:, :, j))) cycle
              l = l + 1                     ! Keep a counter
              temp_elements_uns(l) = ele_B  ! Keep the actual element
!              if ((i .eq. 2) .and. (j .eq. 3) ) print "(i5,a,i5,a,i5,a,i5,a,F19.17,2x,F19.17,8x,F19.17,2x,F19.17,8x,F19.17,2x,F19.17)", mpi_my_id," packing ",l," (ele_B:",ele_B," for process:",j,", element:",ele_val(positionsB, ele_B)
            end do
            allocate(elements_uns(j)%p(l))
            elements_uns(j)%p = temp_elements_uns(1:l)
            deallocate(temp_elements_uns)
! ! ! !            print "(a,i15,a)", "               we are going to SEND  ",l," elements."
            CALL MPI_Isend(l, 1, MPI_INTEGER, j, 0, &
                & MPI_COMM_WORLD, request(k), mpi_my_error)
            allocate(send_buffer(j)%p( l * (positionsB%dim + 1) * 2))
          end if
        end if
      end do
!      write (*,*) ""
!      CALL FLUSH()
!    end if
!  end do

!  CALL FLUSH()
!  CALL MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)

  ! Required GLOBAL communication WAIT
  CALL MPI_Waitall(k, request(1:k), status(:,1:k), mpi_my_error)
!  stop

!  ToDo TODO todo FIX HACK
!  Remove the DO loop and the IF check
!  do i=0,mpi_num_procs-1
!    CALL FLUSH()
!    CALL MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
!    if ( mpi_my_id .eq. i ) then
      do j=0,mpi_num_procs-1
        if (partition_intersection(j) .eqv. .FALSE.) cycle
! ! ! !          print "(i5,a,i2,a,i15,a)", mpi_my_id,", we are RECVing from ",j,":",number_of_elements_to_receive(j)," elements."
          allocate(recv_buffer(j)%p( number_of_elements_to_receive(j) * (positionsB%dim + 1) * 2 ))
          ! JM: Per-processing sizing required
!          allocate(recv_buffer(number_of_elements_to_receive(j)*(positionsB%dim + 1)*2,0:mpi_num_procs-1))
      end do
!    end if
!  end do

!  CALL FLUSH()
!  CALL MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
!  stop

!  allocate(recv_buffer(maxval(parallel_ele_B_array)*(positionsB%dim + 1)*2,0:mpi_num_procs-1))
!  allocate(send_buffer(parallel_ele_B*(positionsB%dim + 1)*2))

!  recv_buffer = -9.0
!  send_buffer = -100.0
  k=1
  l=0
!  do i=0,mpi_num_procs-1
!    CALL FLUSH()
!    CALL MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
!    if ( mpi_my_id .eq. i ) then
      do j=0,mpi_num_procs-1
        k = 1
        if (.NOT. ASSOCIATED(elements_uns(j)%p)) then
!          print "(i5,a,i2,a)", mpi_my_id,", we are NOT sending anything to ",j,". Skipping."
          cycle
        end if
        if ( size(elements_uns(j)%p) .eq. 0 ) then
!          print "(i5,a,i2,a)", mpi_my_id,", we are NOT1 sending anything to ",j,". Skipping."
          cycle
        end if
        allocate(temp_buffer(size(elements_uns(j)%p) * (positionsB%dim + 1) * 2))
!        print "(i5, a,i15,a,i5,a)", mpi_my_id, " we are going to SEND ",size(elements_uns(j)%p)," elements to ",j,"."
        do l=1,size(elements_uns(j)%p)
          tri_B%v = ele_val(positionsB, elements_uns(j)%p(l))
!          if ((i .eq. 2) .and. (j .eq. 3) ) print "(i5,a,i5,a,i5,a,i5,a,F19.17,2x,F19.17,8x,F19.17,2x,F19.17,8x,F19.17,2x,F19.17)", mpi_my_id," packing ",l," (ele_B:",elements_uns(j)%p(l)," for process:",j,", element:",tri_B%v
!!          print "(i5, a,i15,a)", mpi_my_id, " accessing ",l," element."
          temp_buffer(k) = tri_B%v(1,1)
          k = k + 1
          temp_buffer(k) = tri_B%v(2,1)
          k = k + 1
          temp_buffer(k) = tri_B%v(1,2)
          k = k + 1
          temp_buffer(k) = tri_B%v(2,2)
          k = k + 1
          temp_buffer(k) = tri_B%v(1,3)
          k = k + 1
          temp_buffer(k) = tri_B%v(2,3)
          k = k + 1
        end do
        send_buffer(j)%p = temp_buffer
        deallocate(temp_buffer)
!        do ele_B=1,ele_count(positionsB) 
!          if (ele_ownerB(ele_B) /= mpi_my_id) cycle
!          if (elements_uns(l,j) /= ele_B) cycle
!          l = l + 1
!          ! JM: Somewhere in the send loop
!          !if(.not. bboxes_intersect(bbox(ele_val(positionsB, ele_B)), parallel_bbox_a(:, :, proc))) cycle
!    
!          tri_B%v = ele_val(positionsB, ele_B)
!          send_buffer(k) = tri_B%v(1,1)
!          k = k + 1
!          send_buffer(k) = tri_B%v(2,1)
!          k = k + 1
!          send_buffer(k) = tri_B%v(1,2)
!          k = k + 1
!          send_buffer(k) = tri_B%v(2,2)
!          k = k + 1
!          send_buffer(k) = tri_B%v(1,3)
!          k = k + 1
!          send_buffer(k) = tri_B%v(2,3)
!          k = k + 1
!        end do
      end do
!    end if
!  end do

  k = 0
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

      ! Testing, will be removed
!      do j=0,mpi_num_procs-1
!        if (mpi_my_id .eq. j) then
!          recv_buffer(:,j) = 2.0
!          cycle
!        end if
!      end do

      ! Check bounding box intersection
      do j=0,mpi_num_procs-1
        if (mpi_my_id .eq. j) cycle
        if (bboxes_intersect(bbox_a, parallel_bbox_b(:, :, j))) then
!         write (*,*) mpi_my_id," My partition of meshA intersects with ",j," partition of meshB."
!          print "(i5,a,i2,a)", mpi_my_id,", RECVing from ",j,"."
          k = k + 1
! ! !          print "(i5,a,i5,a,i5,a)", mpi_my_id," j:",j,",k:",k,"."
          CALL MPI_Irecv(recv_buffer(j)%p, & 
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
          CALL MPI_Isend(send_buffer(j)%p, &
! ! !             & parallel_ele_B_array(i) * (positionsB%dim + 1), &
! ! !             & parallel_ele_B * (positionsB%dim + 1)*2, &
             & size(send_buffer(j)%p), &
             & MPI_DOUBLE_PRECISION, &
             & j, 0, MPI_COMM_WORLD, request(k), mpi_my_error)
        end if

! ! !    ! For the time being copy Partition B to all Processes
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


!   if ( mpi_my_id .eq. 3 ) then
!     write (*,*) ""
!     tri_A%v = ele_val(positionsA, 1)
!     write (*,*) tri_A%v
! 
!     tri_B%v(1,1) = 0.0000000000000000
!     tri_B%v(2,1) = 0.40000000000166441
!     tri_B%v(1,2) = 7.2188324331636855E-002
!     tri_B%v(2,2) = 0.35587316039890549
!     tri_B%v(1,3) = 0.10967588592745910
!     tri_B%v(2,3) = 0.46046529145073289
!     write (*,*) tri_B%v
!     
!     write (*,*) "BBOX A using ele_val:",bbox(ele_val(positionsA, 1))
!     write (*,*) "BBOX A using %v     :",bbox(tri_A%v)
!     write (*,*) ""
!     write (*,*) "BBOX B using %v     :",bbox(tri_B%v)
!     if(bboxes_partition_intersect(bbox(tri_B%v), bbox(ele_val(positionsA, 1)) )) write (*,*) "bbox(tri_B%v), bbox(ele_val(positionsA, 1))) INTERSECT"
!     if(bboxes_partition_intersect(bbox(ele_val(positionsA, 1)), bbox(tri_B%v)) ) write (*,*) "bbox(ele_val(positionsA, 1)), bbox(tri_B%v)) INTERSECT"
!     if(bboxes_partition_intersect(bbox(tri_B%v), bbox(tri_A%v)))                 write (*,*) "bbox(tri_A%v)               , bbox(tri_B%v)) INTERSECT"
! 
!     call intersect_tris(tri_A, tri_B, trisC, n_trisC)
!     do ele_C=1,n_trisC
!       write (*,*) "area :",triangle_area(trisC(ele_C)%v)
!     end do
! 
! 
! 
!     write (*,*) ""
!     tri_A%v = ele_val(positionsA, 1)
!     write (*,*) tri_A%v
! 
!     tri_B%v(1,1) = 0.0000000000000000
!     tri_B%v(2,1) = 0.50000000000205869
!     tri_B%v(1,2) = 0.00000000000000000
!     tri_B%v(2,2) = 0.40000000000166441
!     tri_B%v(1,3) = 0.10967588592745910
!     tri_B%v(2,3) = 0.46046529145073289
!     write (*,*) tri_B%v
!     
!     write (*,*) "BBOX A using ele_val:",bbox(ele_val(positionsA, 1))
!     write (*,*) "BBOX A using %v     :",bbox(tri_A%v)
!     write (*,*) ""
!     write (*,*) "BBOX B using %v     :",bbox(tri_B%v)
!     if(bboxes_partition_intersect(bbox(tri_B%v), bbox(ele_val(positionsA, 1)) )) write (*,*) "bbox(tri_B%v), bbox(ele_val(positionsA, 1))) INTERSECT"
!     if(bboxes_partition_intersect(bbox(ele_val(positionsA, 1)), bbox(tri_B%v)) ) write (*,*) "bbox(ele_val(positionsA, 1)), bbox(tri_B%v)) INTERSECT"
!     if(bboxes_partition_intersect(bbox(tri_B%v), bbox(tri_A%v)))                 write (*,*) "bbox(tri_A%v)               , bbox(tri_B%v)) INTERSECT"
! 
!     call intersect_tris(tri_A, tri_B, trisC, n_trisC)
!     do ele_C=1,n_trisC
!       write (*,*) "area :",triangle_area(trisC(ele_C)%v)
!     end do
!   end if

!  CALL FLUSH()
!  CALL MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)
!  stop


!  toprint = .TRUE.
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

!        if ( (i .eq. 3) ) then
!        print "(i5,a,F19.17,a)", mpi_my_id,", ORIGINGAL area_parallel (TEMP):",area_parallel,"."
!        end if

      t2 = mpi_wtime()
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
!        if (( ele_A < 5) .and. (mpi_my_id .eq. 0) ) write(*,*) "A: ele_A:",ele_A,", ",tri_A%v
        do ele_B=1,ele_count(positionsB) 
          if(ele_ownerB(ele_B) /= mpi_my_id) cycle
          ! B. Use the libSuperMesh internal triangle intersector (using derived types as input)
          tri_B%v = ele_val(positionsB, ele_B)

          local_iter = local_iter + 1

          call intersect_tris(tri_A, tri_B, trisC, n_trisC)

          do ele_C=1,n_trisC
            area_parallel = area_parallel + triangle_area(trisC(ele_C)%v)
            local_iter_actual = local_iter_actual + 1
          end do
        end do

!        print "(i5,a,i5,a,F19.17,a)", mpi_my_id,", i:",i,", area_parallel (TEMP):",area_parallel,"."
!        CALL FLUSH()

!        do j=0,mpi_num_procs-1
!          if ( ( toprint .eqv. .TRUE. ) .and. (i .eq. 3 ) ) then
!            print "(i5,a,i5,a,i10,a)", mpi_my_id,", Using buffer RECVed from :",j,", RECVed :",number_of_elements_to_receive(j)," elements."
!          end if
!        end do
        toprint = .FALSE.

        do j=0,mpi_num_procs-1
          if ( j .eq. mpi_my_id ) cycle
          if ( partition_intersection(j) .eqv. .FALSE. ) cycle
          if ( number_of_elements_to_receive(j) .eq. 0 ) cycle

          do k=1,number_of_elements_to_receive(j) * (positionsB%dim + 1)*2,6
             tri_B%v(1,1) = recv_buffer(j)%p(k)
             tri_B%v(2,1) = recv_buffer(j)%p(k+1)
             tri_B%v(1,2) = recv_buffer(j)%p(k+2)
             tri_B%v(2,2) = recv_buffer(j)%p(k+3)
             tri_B%v(1,3) = recv_buffer(j)%p(k+4)
             tri_B%v(2,3) = recv_buffer(j)%p(k+5)

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
!        print "(i5,a,i5,a,F19.17,a)", mpi_my_id,", ele_A:",ele_A,", area_parallel:",area_parallel,"."
      end do
!      if ( ele_A < 7 ) write (*,*) ""
      local_time = mpi_wtime() - t0                    ! total
      local_time_intersection_only = mpi_wtime() - t2  ! only intersection
      local_other_time = t2 - t1                       ! other

!      print "(i5,a,e25.17e3)", mpi_my_id, ": intersection time  = ", local_time
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
  write (*,*) ""
  CALL FLUSH()
  CALL MPI_BARRIER(MPI_COMM_WORLD, mpi_my_error)

  ! Gather remote results:
  CALL MPI_Gather(area_parallel, 1, MPI_DOUBLE_PRECISION, &
      & areas_parallel, 1, MPI_DOUBLE_PRECISION, &
      & mpi_my_root, MPI_COMM_WORLD, mpi_my_error)

  CALL MPI_Gather(local_time, 1, MPI_DOUBLE_PRECISION, &
      & times_parallel, 1, MPI_DOUBLE_PRECISION, &
      & mpi_my_root, MPI_COMM_WORLD, mpi_my_error)

  CALL MPI_Gather(local_other_time, 1, MPI_DOUBLE_PRECISION, &
      & other_times_parallel, 1, MPI_DOUBLE_PRECISION, &
      & mpi_my_root, MPI_COMM_WORLD, mpi_my_error)

  CALL MPI_Gather(local_time_intersection_only, 1, MPI_DOUBLE_PRECISION, &
      & times_intersection_only_parallel, 1, MPI_DOUBLE_PRECISION, &
      & mpi_my_root, MPI_COMM_WORLD, mpi_my_error)

  CALL MPI_Gather(local_mesh, 1, MPI_DOUBLE_PRECISION, &
      & mesh_parallel, 1, MPI_DOUBLE_PRECISION, &
      & mpi_my_root, MPI_COMM_WORLD, mpi_my_error)

  CALL MPI_Gather(local_iter, 1, MPI_INTEGER, &
      & iters_parallel, 1, MPI_INTEGER, &
      & mpi_my_root, MPI_COMM_WORLD, mpi_my_error)

  CALL MPI_Gather(local_iter_actual, 1, MPI_INTEGER, &
      & iter_actual_parallel, 1, MPI_INTEGER, &
      & mpi_my_root, MPI_COMM_WORLD, mpi_my_error)

  CALL MPI_REDUCE(parallel_ele_A, local_sum_a, 1, MPI_INT, MPI_SUM, mpi_my_root, MPI_COMM_WORLD, mpi_my_error)

  CALL MPI_REDUCE(parallel_ele_B, local_sum_b, 1, MPI_INT, MPI_SUM, mpi_my_root, MPI_COMM_WORLD, mpi_my_error)

  if ( mpi_my_id .eq. 0 ) then 
!    write (*,*) times_parallel
!    write (*,*) areas_parallel
!    print "(i5,a,i15,a,i15,a)", mpi_my_id, ": Ele A:",serial_ele_A,", Ele B:",serial_ele_B,"."
    print "(i5,a,F19.15,a)", mpi_my_id, ": Total serial read time    : ", mesh_serial," ."
    print "(i5,a,F19.15,a)", mpi_my_id, ": Total parallel read time  : ", sum(mesh_parallel)," ."
    print "(i5,a,F19.15,a)", mpi_my_id, ": Max parallel read time    : ", maxval(mesh_parallel)," ."
    write(*,*) ""
    print "(i5,a,F19.15,a)", mpi_my_id, ": Total serial intersection time  : ", time_serial," ."
    print "(i5,a,F19.15,a)", mpi_my_id, ": Total parallel intersection time: ", sum(times_intersection_only_parallel)," ."
    print "(i5,a,F19.15,a)", mpi_my_id, ": Max parallel intersection time  : ", maxval(times_intersection_only_parallel)," ."
    write(*,*) ""
    print "(i5,a,F19.15,a)", mpi_my_id, ": Total serial time               : ", time_serial," ."
    print "(i5,a,F19.15,a)", mpi_my_id, ": Total parallel time             : ", sum(times_parallel)," ."
    print "(i5,a,F19.15,a)", mpi_my_id, ": Max parallel time               : ", maxval(times_parallel)," ."
    write(*,*) ""
    print "(i5,a,F19.15,a)", mpi_my_id, ": Total other parallel time       : ", sum(other_times_parallel)," ."
    print "(i5,a,F19.15,a)", mpi_my_id, ": Max other parallel time         : ", maxval(other_times_parallel)," ."
    write(*,*) ""
    print "(i5,a,F19.15,a)", mpi_my_id, ": Total serial intersection area     : ", area_serial," ."
    print "(i5,a,F19.15,a)", mpi_my_id, ": Total parallel intersection area   : ", sum(areas_parallel)," ."
!    print "(i5,a,i15,a,i15,a)", mpi_my_id, ": Serial iters         :",serial_local_iter,"."
!    print "(i5,a,i15,a,i15,a)", mpi_my_id, ": Parallel iters       :",sum(iters_parallel),"."
!    print "(i5,a,i15,a,i15,a)", mpi_my_id, ": Serial actual iters  :",serial_local_iter_actual,"."
!    print "(i5,a,i15,a,i15,a)", mpi_my_id, ": Parallel actual iters:",sum(iter_actual_parallel),"."

!    write (*,*) ""
!    write (*,*) "local_iter_actual_parallel:",iter_actual_parallel

    fail = fnequals(sum(areas_parallel), area_serial, tol = tol)
    call report_test("[test_parallel_partition_ab areas]", fail, .FALSE., "Should give the same areas of intersection")
    if (fail) then
      print "(a,e25.17e3,a,e25.17e3,a)", ": Total parallel intersection area:", sum(areas_parallel),&
                                       & ", total serial intersection area:",area_serial,"."
      write (*,*) "areas_parallel:",areas_parallel
      write (*,*) "iters_parallel:",iters_parallel
      print "(a,i15,a,i15,a)", "Serial iters:",serial_local_iter,", parallel iters:",sum(iters_parallel),"."
    end if

    fail = ( sum(iter_actual_parallel) .ne. serial_local_iter_actual)
    call report_test("[test_parallel_partition_ab iterations]", fail, .FALSE., "Should give the same number of iterations")
    if (fail) then
      print "(a,i15,a,i15,a)", ": Total parallel actual iterations:", sum(iter_actual_parallel),&
                            & ", total serial actual iterations  :",serial_local_iter_actual,"."
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

  do i=0,size(elements_uns(:))-1
    if ( .NOT. ASSOCIATED(elements_uns(i)%p) ) cycle
    deallocate(elements_uns(i)%p)
  end do
  deallocate(elements_uns)

  do i=0,size(send_buffer(:))-1
    if ( .NOT. ASSOCIATED(send_buffer(i)%p) ) cycle
    deallocate(send_buffer(i)%p)
  end do
  deallocate(send_buffer)

  do i=0,size(recv_buffer(:))-1
    if ( .NOT. ASSOCIATED(recv_buffer(i)%p) ) cycle
    deallocate(recv_buffer(i)%p)
  end do
  deallocate(recv_buffer)

  deallocate(areas_parallel, times_parallel)
  deallocate(parallel_bbox_a, parallel_bbox_b)
  deallocate(parallel_ele_B_array)
  deallocate(request, status, partition_intersection, times_intersection_only_parallel)

  call cintersection_finder_reset(nnodes)

end subroutine test_parallel_partition_ab
