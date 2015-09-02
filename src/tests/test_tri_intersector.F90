subroutine test_tri_intersector

  use libsupermesh_construction
  use libsupermesh_fields
  use libsupermesh_read_triangle_2
  use libsupermesh_tri_intersection_module
  use libsupermesh_unittest_tools
  
  implicit none

  type(vector_field) :: positionsA, positionsB
  integer :: ele_A, ele_B, ele_C
  real :: area_D_libwm, area_A_libwm_intersect, area_B_fort, area_C_fort_public, &
     &  area_E_intersect_elements, area_F_intersect_elements, area_G_intersect_elements, &
     & area_H_intersect_old
  logical :: fail, libWM_tri_failed
  
  type(tri_type) :: triA, triB
  type(tri_type), dimension(tri_buf_size) :: trisC
  integer :: ntests, n_trisC, i, nonods, totele, n_trisC_pre
  integer, dimension(:, :), allocatable :: ndglno
  real, dimension(2, tri_buf_size) :: nodesC
  integer, dimension(3, tri_buf_size) :: ndglnoC
  real, dimension(2, 3, tri_buf_size) ::  trisC_real
  type(vector_field) :: intersection
  type(mesh_type) :: intersection_mesh, new_mesh

  integer, parameter :: dim = 2, loc = 3
  
  positionsA = read_triangle_files("data/plcA", dim = dim)
  positionsB = read_triangle_files("data/plcB", dim = dim)

  call cintersector_set_dimension(dim)

  do ele_A=1,ele_count(positionsA)
    do ele_B=1,ele_count(positionsB)

      triA%v = ele_val(positionsA, ele_A)
      triB%v = ele_val(positionsB, ele_B)

      ! A. Use libWM without creating temporary vector fields.
      call intersect_tris(triA%v, triB%v, nodesC, ndglnoC, n_trisC)
      area_A_libwm_intersect = 0.0
      do ele_C=1,n_trisC
        area_A_libwm_intersect = area_A_libwm_intersect + triangle_area(nodesC(:, ndglnoC(:, ele_C)))
      end do
      n_trisC_pre = n_trisC

      ! B. Use the libSuperMesh internal triangle intersector (using derived types as input)
      call intersect_tris(triA, triB, trisC, n_trisC)
      area_B_fort = 0.0
      do ele_C=1,n_trisC
        area_B_fort = area_B_fort + triangle_area(trisC(ele_C)%v)
      end do
      
!       if (n_trisC_pre .ne. n_trisC) then
!         write(*,*) "LibWM returned:",n_trisC_pre,", and tris returned:",n_trisC,"."
!         do ele_C=1,n_trisC_pre
!           write(*,*) "libWM: area:",triangle_area(nodesC(:, ndglnoC(:, ele_C))),"."
!           temp_area_A_libwm_intersect = temp_area_A_libwm_intersect + triangle_area(nodesC(:, ndglnoC(:, ele_C)))
!         end do
!         do ele_C=1,n_trisC
!           write(*,*) "tris : area:",triangle_area(trisC(ele_C)%v),"."
!           temp_area_B_fort = temp_area_B_fort + triangle_area(trisC(ele_C)%v)
!         end do
! 
!         if ( abs(temp_area_B_fort - temp_area_A_libwm_intersect) > epsilon(0.0) ) then
!           write(*,*) "Input triangles"
!           write(*,*) "A: ele_A:",ele_A,", triA:",ele_val(positionsA, ele_A),"."
!           write(*,*) "B: ele_B:",ele_B,", triB:",ele_val(positionsB, ele_B),"."
!           write(*,*) "Results:"
!           do ele_C=1,n_trisC_pre
!             write(*,*) "libWM: area:",triangle_area(nodesC(:, ndglnoC(:, ele_C))),", nodes:",nodesC(:, ndglnoC(:, ele_C)),"."
!           end do
!           do ele_C=1,n_trisC
!             write(*,*) "tris : area:",triangle_area(trisC(ele_C)%v),",nodes:",trisC(ele_C)%v,"."
!           end do
!           call exit(1)
!         end if
!         libWM_tri_failed = .TRUE.
!       end if

      ! C. Use the libSuperMesh internal triangle intersector (using only reals as input)
      call intersect_tris(triA%v, triB%v, trisC_real, n_trisC)
      area_C_fort_public = 0.0
      do ele_C=1,n_trisC
        area_C_fort_public = area_C_fort_public + triangle_area(trisC_real(:,:,ele_C))
      end do

      ! D. Use libWM directly and create temporary vector field
      call cintersector_set_input(ele_val(positionsA, ele_A), triB%v, dim, loc)
      call cintersector_drive
      call cintersector_query(nonods, totele)
      call allocate(intersection_mesh, dim, nonods, totele, loc)
      intersection_mesh%continuity = -1
      call allocate(intersection, dim, intersection_mesh)
      if (nonods > 0) then
        allocate(ndglno(loc, totele))
        call cintersector_get_output(nonods, totele, dim, loc, intersection%val, ndglno)
        intersection_mesh%ndglno = reshape(ndglno, (/loc * totele/))
        deallocate(ndglno)
      end if
      area_D_libwm = 0.0
      do ele_C=1,totele
        area_D_libwm = area_D_libwm + triangle_area(ele_val(intersection, ele_C))
      end do
      call deallocate(intersection_mesh)
      call deallocate(intersection)

      ! E. Use libWM directly and create temporary vector field
      intersection = intersect_elements_old(triA%v, &
        triB%v, loc, dim)
      area_E_intersect_elements = 0.0
      do ele_C=1,ele_count(intersection)
        area_E_intersect_elements = area_E_intersect_elements + triangle_area(ele_val(intersection, ele_C))
      end do
      call deallocate(intersection)

      ! F. Use the new intersect_elements and do NOT create vector field
      call intersect_elements(triA%v, triB%v, n_trisC, trisC_real)
      area_F_intersect_elements = 0.0
      do ele_C=1,n_trisC
        area_F_intersect_elements = area_F_intersect_elements + triangle_area(trisC_real(:,:,ele_C))
      end do

      ! G. Use the new intersect_elements and DO create vector field
      call intersect_elements(triA%v, triB%v, n_trisC, trisC_real)
      call allocate(new_mesh, dim, n_trisC * loc, n_trisC, loc)

      if ( n_trisC > 0 ) then
        new_mesh%ndglno = (/ (i, i=1,loc * n_trisC) /)
        new_mesh%continuity = -1
      end if

      call allocate(intersection, dim, new_mesh)
      if ( n_trisC > 0 ) then
        do i = 1, n_trisC
          call set(intersection, ele_nodes(intersection, i), trisC_real(:,:,i))
        end do
      end if
      call deallocate(new_mesh)
      area_G_intersect_elements = 0.0
      do ele_C=1,ele_count(intersection)
        area_G_intersect_elements = area_G_intersect_elements + triangle_area(ele_val(intersection, ele_C))
      end do
      call deallocate(intersection)
      
      ! H. Use the *OLD* libSuperMesh internal triangle intersector (using derived types as input)
      call intersect_tris_dt_old(triA, triB, trisC, n_trisC)
      do ele_C=1,n_trisC
          area_H_intersect_old = area_H_intersect_old + triangle_area(trisC(ele_C)%v)
      end do

!       if (libWM_tri_failed .eqv. .TRUE.) then
!         write(*,*) "LibWM returned:",n_trisC_pre,", and tris_old returned:",n_trisC,"."
!         do ele_C=1,n_trisC_pre
!           write(*,*) "libWM  : area:",triangle_area(nodesC(:, ndglnoC(:, ele_C))),"."
!           temp_area_A_libwm_intersect = temp_area_A_libwm_intersect + triangle_area(nodesC(:, ndglnoC(:, ele_C)))
!         end do
!         do ele_C=1,n_trisC
!           write(*,*) "tris_old: area:",triangle_area(trisC(ele_C)%v),"."
!           temp_area_H_intersect_old = temp_area_H_intersect_old + triangle_area(trisC(ele_C)%v)
!         end do
! 
!         if ( abs(temp_area_H_intersect_old - temp_area_A_libwm_intersect) > epsilon(0.0) ) then
!           write(*,*) "Input triangles"
!           write(*,*) "A: ele_A:",ele_A,", triA:",ele_val(positionsA, ele_A),"."
!           write(*,*) "B: ele_B:",ele_B,", triB:",ele_val(positionsB, ele_B),"."
!           write(*,*) "Results:"
!           do ele_C=1,n_trisC_pre
!             write(*,*) "libWM:    area:",triangle_area(nodesC(:, ndglnoC(:, ele_C))),", nodes:",nodesC(:, ndglnoC(:, ele_C)),"."
!           end do
!           do ele_C=1,n_trisC
!             write(*,*) "tris_old: area:",triangle_area(trisC(ele_C)%v),",nodes:",trisC(ele_C)%v,"."
!           end do
!           
!           if (n_trisC_pre .ne. n_trisC) then
! !            call exit(1)
!           end if
!           write (*,*)
!         end if
!       end if

      fail = (area_A_libwm_intersect .fne. area_B_fort) &
         .OR. (area_E_intersect_elements .fne. area_C_fort_public ) &
         .OR. (area_E_intersect_elements .fne. area_B_fort) &
         .OR. (area_A_libwm_intersect .fne. area_D_libwm) &
         .OR. (area_A_libwm_intersect .fne. area_C_fort_public ) &
         .OR. (area_F_intersect_elements .fne. area_A_libwm_intersect ) &
         .OR. (area_F_intersect_elements .fne. area_D_libwm ) &
         .OR. (area_G_intersect_elements .fne. area_B_fort) &
         .OR. (area_G_intersect_elements .fne. area_D_libwm )
      call report_test("[tri_intersector areas]",fail, .false., "Should give the same areas of intersection")
      if ( fail .eqv. .TRUE. ) then
        write (*,*) "[tri_intersector areas] area_A_libwm_intersect:",area_A_libwm_intersect, &
          ", area_B_fort:", area_B_fort,", area_C_fort_public:",area_C_fort_public,&
          ", area_D_libwm:",area_D_libwm, &
          ", area_E_intersect_elements:",area_E_intersect_elements, &
          ", area_F_intersect_elements:",area_F_intersect_elements, &
          ", area_G_intersect_elements:",area_G_intersect_elements,"."
      end if

      libWM_tri_failed = .FALSE.
    end do
  end do

  call cintersection_finder_reset(ntests)
  
  call deallocate(positionsA)
  call deallocate(positionsB)

end subroutine test_tri_intersector
