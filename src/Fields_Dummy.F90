module fields_dummy

  use libsupermesh_global_parameters, only : FIELD_NAME_LEN
  use libsupermesh_integer_set_module

  implicit none

  private

  public :: mesh_type, vector_field, eelist_type, allocate, deallocate, extract_eelist, ele_count, node_count, ele_val, ele_loc, node_val, row_m_ptr

  type mesh_type
    character(len = FIELD_NAME_LEN) :: name
    integer :: nnodes
    integer :: nelements
    integer :: loc
    integer, dimension(:), pointer :: ndglno
    integer :: continuity
    type(eelist_type), pointer :: eelist
  end type mesh_type

  type vector_field
    character(len = FIELD_NAME_LEN) :: name
    type(mesh_type), pointer :: mesh
    integer :: dim
    real, dimension(:, :), pointer :: val
  end type vector_field

  type eelist_type
    integer, dimension(:, :), pointer :: v
    integer, dimension(:), pointer :: n
  end type eelist_type

  interface allocate
    module procedure allocate_mesh, allocate_vector_field
  end interface allocate

  interface deallocate
    module procedure deallocate_mesh, deallocate_vector_field, deallocate_eelist
  end interface deallocate

  interface extract_eelist
    module procedure extract_eelist_mesh, extract_eelist_vector_field
  end interface extract_eelist

  interface ele_count
    module procedure ele_count_vector_field
  end interface ele_count

  interface node_count
    module procedure node_count_vector_field
  end interface node_count

  interface ele_loc
    module procedure ele_loc_vector_field
  end interface ele_loc

contains

  subroutine allocate_mesh(mesh, nnodes, nelements, loc, name)
    type(mesh_type), intent(out) :: mesh
    integer, intent(in) :: nnodes
    integer, intent(in) :: nelements
    integer, intent(in) :: loc
    character(len = *), optional, intent(in) :: name

    mesh%nnodes = nnodes
    mesh%nelements = nelements
    mesh%loc = loc
    allocate(mesh%ndglno(nelements * loc))
    if(present(name)) then
      mesh%name = name
    else
      mesh%name = "mesh"
    end if
    nullify(mesh%eelist)
    continuity = 0
    
  end subroutine allocate_mesh

  subroutine deallocate_mesh(mesh)
    type(mesh_type), intent(inout) :: mesh

    deallocate(mesh%ndglno)
    if(associated(mesh%eelist)) then
      call deallocate(mesh%eelist)
      deallocate(mesh%eelist)
    end if

  end subroutine deallocate_mesh

  subroutine allocate_vector_field(field, dim, mesh, name)
    type(vector_field), intent(out) :: field
    integer, intent(in) :: dim
    type(mesh_type), target, intent(in) :: mesh
    character(len = *), optional, intent(in) :: name

    field%dim = dim
    allocate(field%val(dim, mesh%nnodes))
    field%mesh => mesh
    if(present(name)) then
      field%name = name
    else
      field%name = "field"
    end if
    
  end subroutine allocate_vector_field

  subroutine deallocate_vector_field(field)
    type(vector_field), intent(inout) :: field

    deallocate(field%val)

  end subroutine deallocate_vector_field

  subroutine mesh_eelist(nnodes, enlist, eelist)
    integer, intent(in) :: nnodes
    ! loc x nelements
    integer, dimension(:, :), intent(in) :: enlist
    type(eelist_type), intent(out) :: eelist

    integer :: ele, i, j, k, loc, nelements, nneigh_ele
    type(integer_set), dimension(:), allocatable :: nelist

    integer, dimension(:), allocatable :: seen, seen_eles
    integer :: nseen_eles

    type iarr
      integer, dimension(:), pointer :: v
      integer :: n
    end type iarr
    type(iarr), dimension(:), allocatable :: inelist

    loc = size(enlist, 1)
    nelements = size(enlist, 2)

    ! Construct the node element list
    allocate(nelist(nnodes))
    call allocate(nelist)
    do i = 1, nelements
      do j = 1, loc
        call insert(nelist(enlist(j, i)), i)
      end do
    end do
    ! Copy the node element list into a ragged array for more efficient access
    allocate(inelist(nnodes))
    do i = 1, nnodes
      inelist(i)%n = key_count(nelist(i))
      allocate(inelist(i)%v(inelist(i)%n))
      do j = 1, inelist(i)%n
        inelist(i)%v(j) = fetch(nelist(i), j)
      end do
      call deallocate(nelist(i))
    end do
    deallocate(nelist)

    allocate(seen_eles(nelements), seen(nelements))
    seen = 0
    nneigh_ele = loc  ! Linear simplex assumption here (works for more general
                      ! cases, but less efficient)
    allocate(eelist%v(nneigh_ele, nelements), eelist%n(nelements))
    eelist%n = 0
    do i = 1, nelements
      nseen_eles = 0
      loc_loop: do j = 1, loc
        do k = 1, inelist(enlist(j, i))%n
          ele = inelist(enlist(j, i))%v(k)
          if(ele /= i) then
            seen(ele) = seen(ele) + 1
            if(seen(ele) == 1) then
              nseen_eles = nseen_eles + 1
              seen_eles(nseen_eles) = ele
            else if(seen(ele) == 2) then
              eelist%n(i) = eelist%n(i) + 1
#ifdef NDEBUG
              eelist%v(eelist%n(i), i) = ele
              if(eelist%n(i) == nneigh_ele) exit loc_loop
#else
              if(eelist%n(i) > nneigh_ele) then
                write(0, "(a)") "Invalid connectivity"
                stop 1
              end if
              eelist%v(eelist%n(i), i) = ele
            else
              write(0, "(a)") "Invalid connectivity"
              stop 1
#endif
            end if
          end if
        end do
      end do loc_loop
      do j = 1, nseen_eles
        seen(seen_eles(j)) = 0
      end do
    end do
    deallocate(seen_eles, seen)

    do i = 1, nnodes
      deallocate(inelist(i)%v)
    end do
    deallocate(inelist)

  end subroutine mesh_eelist

  subroutine deallocate_eelist(eelist)
    type(eelist_type), intent(inout) :: eelist

    deallocate(eelist%v, eelist%n)

  end subroutine deallocate_eelist

  function extract_eelist_mesh(mesh) result(eelist)
    type(mesh_type), intent(inout) :: mesh

    type(eelist_type), pointer :: eelist

    if(.not. associated(mesh%eelist)) then
      allocate(mesh%eelist)
      call mesh_eelist(mesh%nnodes, reshape(mesh%ndglno, (/mesh%loc, mesh%nelements/)), mesh%eelist)
    end if
    eelist => mesh%eelist

  end function extract_eelist_mesh

  function extract_eelist_vector_field(field) result(eelist)
    type(vector_field), intent(inout) :: field

    type(eelist_type), pointer :: eelist

    eelist => extract_eelist(field%mesh)

  end function extract_eelist_vector_field

  pure function ele_count_vector_field(field) result(ele_count)
    type(vector_field), intent(in) :: field

    integer :: ele_count

    ele_count = field%mesh%nelements

  end function ele_count_vector_field

  pure function node_count_vector_field(field) result(node_count)
    type(vector_field), intent(in) :: field

    integer :: node_count

    node_count = field%mesh%nnodes

  end function node_count_vector_field

  pure function ele_val(field, ele)
    type(vector_field), intent(in) :: field
    integer, intent(in) :: ele

    real, dimension(field%dim, field%mesh%loc) :: ele_val

    ele_val = field%val(:, field%mesh%ndglno((ele - 1) * field%mesh%loc + 1:ele * field%mesh%loc))

  end function ele_val

  pure function ele_loc_vector_field(field, ele) result(ele_loc)
    type(vector_field), intent(in) :: field
    integer, intent(in) :: ele

    integer :: ele_loc

    ele_loc = field%mesh%loc

  end function ele_loc_vector_field

  pure function node_val(field, node)
    type(vector_field), intent(in) :: field
    integer, intent(in) :: node

    real, dimension(field%dim) :: node_val

    node_val = field%val(:, node)

  end function node_val

  function row_m_ptr(eelist, i)
    type(eelist_type), intent(in) :: eelist
    integer, intent(in) :: i

    integer, dimension(:), pointer :: row_m_ptr

    row_m_ptr => eelist%v(:eelist%n(i), i)

  end function row_m_ptr

end module fields_dummy