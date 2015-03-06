#include "fdebug.h"

module libsupermesh_quicksort

  use libsupermesh_fldebug
  use libsupermesh_global_parameters, only : real_4

  implicit none

  private
  
!  public :: qsort, sort, count_unique, inverse_permutation, apply_permutation, apply_reverse_permutation
  public :: count_unique

contains

  function count_unique(int_array) result(unique)
    !!< Count the unique entries in the supplied array of integers

    integer, dimension(:), intent(in) :: int_array

    integer :: unique

    integer :: i
    integer, dimension(size(int_array)) :: permutation

    call qsort(int_array, permutation)
    
    unique = 0
    if(size(int_array) > 0) then
      unique = unique + 1
    end if
    do i = 2, size(int_array)
      if(int_array(permutation(i)) == int_array(permutation(i - 1))) cycle
      unique = unique + 1
    end do

  end function count_unique
  
end module libsupermesh_quicksort
