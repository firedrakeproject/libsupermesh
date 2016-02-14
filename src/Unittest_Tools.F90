#include "libsupermesh_configuration.h"

module libsupermesh_unittest_tools
  ! Utility functions for unit testing

  use libsupermesh_debug_parameters, only : debug_log_unit

  implicit none

  private
  
  public :: report_test, operator(.feq.), operator(.fne.), fequals, fnequals

  interface operator(.feq.)
    module procedure fequals_scalar_op, fequals_array_op, &
      & fequals_array_scalar_op, fequals_matrix_op, fequals_matrix_scalar_op
  end interface operator(.feq.)

  interface operator(.fne.)
    module procedure fnequals_scalar_op, fnequals_array_op, &
      & fnequals_array_scalar_op, fnequals_matrix_op, fnequals_matrix_scalar_op
  end interface operator(.fne.)

  interface fequals
    module procedure fequals_scalar, fequals_array, fequals_array_scalar, &
      & fequals_matrix, fequals_matrix_scalar
  end interface fequals

  interface fnequals
    module procedure fnequals_scalar, fnequals_array, fnequals_array_scalar, &
      & fnequals_matrix, fnequals_matrix_scalar 
  end interface fnequals

contains

  subroutine report_test(title, fail, warn, msg)
    !!< This is the subroutine used by unit tests to report the output of
    !!< a test case.

    !! title: the name of the test case.
    character(len = *), intent(in) :: title
    !! msg: an explanatory message printed if the test case fails.
    character(len = *), intent(in) :: msg
    !! Has the test case failed, or triggered a warning? Set fail or warn to .true. if so.
    logical, intent(in) :: fail, warn

    if(fail) then
      write(debug_log_unit, "(a)") "(Fail: " // trim(title) // "; error: " // trim(msg) // ")"
    else if(warn) then
      write(debug_log_unit, "(a)") "(Warn: " // trim(title) // "; error: " // trim(msg) // ")"
    else
      write(debug_log_unit, "(a)") "(Pass: " // trim(title) // ")"
    end if

  end subroutine report_test

  pure function fequals_scalar_op(float1, float2) result(equals)
    real, intent(in) :: float1
    real, intent(in) :: float2

    logical :: equals

    equals = fequals(float1, float2)

  end function fequals_scalar_op

  pure function fnequals_scalar_op(float1, float2) result(nequals)
    real, intent(in) :: float1
    real, intent(in) :: float2

    logical :: nequals

    nequals = fnequals(float1, float2)

  end function fnequals_scalar_op

  pure function fequals_array_op(array1, array2) result(equals)
    real, dimension(:), intent(in) :: array1
    real, dimension(size(array1)), intent(in) :: array2

    logical :: equals

    equals = fequals(array1, array2)

  end function fequals_array_op

  pure function fnequals_array_op(array1, array2) result(nequals)
    real, intent(in), dimension(:) :: array1, array2

    logical :: nequals

    nequals = fnequals(array1, array2)

  end function fnequals_array_op

  pure function fequals_array_scalar_op(array1, float2) result(equals)
    real, dimension(:), intent(in) :: array1
    real, intent(in) :: float2

    logical :: equals

    equals = fequals(array1, float2)

  end function fequals_array_scalar_op

  pure function fnequals_array_scalar_op(array1, float2) result(nequals)
    real, dimension(:), intent(in) :: array1
    real, intent(in) :: float2

    logical :: nequals

    nequals = fnequals(array1, float2)

  end function fnequals_array_scalar_op
  
  pure function fequals_matrix_op(mat1, mat2) result(equals)
    real, dimension(:, :), intent(in) :: mat1
    real, dimension(size(mat1, 1), size(mat1, 2)), intent(in) :: mat2
    
    logical :: equals

    equals = fequals(mat1, mat2)

  end function fequals_matrix_op
  
  pure function fnequals_matrix_op(mat1, mat2) result(nequals)
    real, dimension(:, :), intent(in) :: mat1
    real, dimension(size(mat1, 1), size(mat1, 2)), intent(in) :: mat2
    
    logical :: nequals

    nequals = fnequals(mat1, mat2)

  end function fnequals_matrix_op
  
  pure function fequals_matrix_scalar_op(mat1, float2) result(equals)
    real, dimension(:, :), intent(in) :: mat1
    real, intent(in) :: float2
    
    logical :: equals

    equals = fequals(mat1, float2)

  end function fequals_matrix_scalar_op
  
  pure function fnequals_matrix_scalar_op(mat1, float2) result(nequals)
    real, dimension(:, :), intent(in) :: mat1
    real, intent(in) :: float2
    
    logical :: nequals

    nequals = fnequals(mat1, float2)

  end function fnequals_matrix_scalar_op

  pure function fequals_scalar(float1, float2, tol) result(equals)
    real, intent(in) :: float1
    real, intent(in) :: float2
    real, intent(in), optional :: tol
    
    logical :: equals

    real :: eps
    
    if(present(tol)) then
      eps = abs(tol)
    else
      eps = 1.0D2 * max(epsilon(0.0D0), spacing(float1), spacing(float2))
    end if
    equals = abs(float1 - float2) < eps

  end function fequals_scalar

  pure function fnequals_scalar(float1, float2, tol) result(nequals)
    real, intent(in) :: float1
    real, intent(in) :: float2
    real, optional, intent(in) :: tol

    logical :: nequals

    nequals = .not. fequals(float1, float2, tol = tol)

  end function fnequals_scalar

  pure function fequals_array(array1, array2, tol) result(equals)
    real, dimension(:), intent(in) :: array1
    real, dimension(size(array1)), intent(in) :: array2
    real, intent(in), optional :: tol

    logical :: equals

    integer :: i

    do i = 1, size(array1)
      if(fnequals(array1(i), array2(i), tol = tol)) then
        equals = .false.
        return
      end if
    end do
    equals = .true.

  end function fequals_array

  pure function fnequals_array(array1, array2, tol) result(nequals)
    real, dimension(:), intent(in) :: array1
    real, dimension(size(array1)), intent(in) :: array2
    real, intent(in), optional :: tol

    logical :: nequals

    nequals = .not. fequals(array1, array2, tol = tol)

  end function fnequals_array

  pure function fequals_array_scalar(array1, float2, tol) result(equals)
    real, dimension(:), intent(in) :: array1
    real, intent(in) :: float2
    real, intent(in), optional :: tol

    logical :: equals
    
    integer :: i

    do i = 1, size(array1)
      if(fnequals(array1(i), float2, tol = tol)) then
        equals = .false.
        return
      end if
    end do
    equals = .true.

  end function fequals_array_scalar

  pure function fnequals_array_scalar(array1, float2, tol) result(nequals)
    real, dimension(:), intent(in) :: array1
    real, intent(in) :: float2
    real, intent(in), optional :: tol

    logical :: nequals

    nequals = .not. fequals(array1, float2, tol = tol)

  end function fnequals_array_scalar

  pure function fequals_matrix(mat1, mat2, tol) result(equals)
    real, dimension(:, :), intent(in) :: mat1
    real, dimension(size(mat1, 1), size(mat1, 2)), intent(in) :: mat2
    real, optional, intent(in) :: tol
    
    logical :: equals

    integer :: i

    do i = 1, size(mat1, 1)
      if(fnequals(mat1(i, :), mat2(i, :), tol = tol)) then
        equals = .false.
        return
      end if
    end do
    equals = .true.

  end function fequals_matrix
  
  pure function fnequals_matrix(mat1, mat2, tol) result(nequals)
    real, dimension(:, :), intent(in) :: mat1
    real, dimension(size(mat1, 1), size(mat1, 2)), intent(in) :: mat2
    real, optional, intent(in) :: tol
    
    logical :: nequals

    nequals = .not. fequals(mat1, mat2, tol = tol)

  end function fnequals_matrix

  pure function fequals_matrix_scalar(mat1, float2, tol) result(equals)
    real, dimension(:, :), intent(in) :: mat1
    real, intent(in) :: float2
    real, optional, intent(in) :: tol
    
    logical :: equals

    integer :: i

    do i = 1, size(mat1, 1)
      if(fnequals(mat1(i, :), float2, tol = tol)) then
        equals = .false.
        return
      end if
    end do
    equals = .true.

  end function fequals_matrix_scalar
  
  pure function fnequals_matrix_scalar(mat1, float2, tol) result(nequals)
    real, dimension(:, :), intent(in) :: mat1
    real, intent(in) :: float2
    real, optional, intent(in) :: tol
    
    logical :: nequals

    nequals = .not. fequals(mat1, float2, tol = tol)

  end function fnequals_matrix_scalar

end module libsupermesh_unittest_tools
