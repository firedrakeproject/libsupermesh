module libsupermesh_parallel_supermesh

  implicit none

contains

  subroutine parallel_supermesh(positions_a, enlist_a, positions_b, enlist_b, donor_ele_data, unpack_donor_ele_data, target_ele_data, unpack_target_ele_data, intersection_calculation)
    ! dim x nnodes_a
    real, dimension(:, :), intent(in) :: positions_a
    ! loc_a x nelements_a
    integer, dimension(:, :), intent(in) :: enlist_a
    ! dim x nnodes_b
    real, dimension(:, :), intent(in) :: positions_b
    ! loc_b x nelements_b
    integer, dimension(:, :), intent(in) :: enlist_b
    interface
      subroutine donor_ele_data(eles, data, ndata)
        use iso_c_binding, only : c_ptr
        implicit none
        integer, dimension(:), intent(in) :: eles
        type(c_ptr), intent(out) :: data
        integer, intent(out) :: ndata
      end subroutine donor_ele_data
      
      subroutine unpack_donor_ele_data(ele, data, ndata, ele_data, nele_data)
        use iso_c_binding, only : c_ptr
        implicit none
        integer, intent(in) :: ele
        type(c_ptr), intent(in) :: data
        integer, intent(in) :: ndata
        type(c_ptr), intent(out) :: ele_data
        integer, intent(out) :: nele_data
      end subroutine unpack_donor_ele_data
      
      subroutine target_ele_data(eles, data, ndata)
        use iso_c_binding, only : c_ptr
        implicit none
        integer, dimension(:), intent(in) :: eles
        type(c_ptr), intent(out) :: data
        integer, intent(out) :: ndata
      end subroutine target_ele_data

      subroutine unpack_target_ele_data(ele, data, ndata, ele_data, nele_data)
        use iso_c_binding, only : c_ptr
        implicit none
        integer, intent(in) :: ele
        type(c_ptr), intent(in) :: data
        integer, intent(in) :: ndata
        type(c_ptr), intent(out) :: ele_data
        integer, intent(out) :: nele_data
      end subroutine unpack_target_ele_data

      subroutine intersection_calculation(positions_c, ele_a, proc_a, ele_b, ele_data_a, ele_ndata_a, ele_data_b, ele_ndata_b)
        use iso_c_binding, only : c_ptr
        implicit none
        ! dim x (loc_c x nelements_c)
        real, dimension(:, :), intent(in) :: positions_c
        integer, intent(in) :: ele_a
        integer, intent(in) :: proc_a
        integer, intent(in) :: ele_b
        type(c_ptr), intent(in) :: ele_data_a
        integer, intent(in) :: ele_ndata_a
        type(c_ptr), intent(in) :: ele_data_b
        integer, intent(in) :: ele_ndata_b
      end subroutine intersection_calculation
    end interface

    ! 1. Communicate bounding box data for donor and target

    ! 2. Use bounding box data to cull donor mesh

    ! 3. Pack non-culled mesh data for communication, calling user specified data functions

    ! 4. Communicate mesh and mesh data

    ! 5. Supermesh and call user specified element unpack and calculation functions
    
  end subroutine parallel_supermesh

end module libsupermesh_parallel_supermesh