#include "libsupermesh_debug.h"

module libsupermesh_interval_intersection

  implicit none
  
  private
  
  public :: intersect_intervals
  
  integer, parameter, public :: interval_buf_size = 1
  
  interface intersect_intervals
    module procedure intersect_intervals_rank_1, intersect_intervals_rank_2, &
      & intersect_intervals_rank_3
  end interface intersect_intervals
  
contains

  pure subroutine intersect_intervals_rank_1(interval_a, interval_b, intervals_c, n_intervals_c)
    real, dimension(2), intent(in) :: interval_a
    real, dimension(2), intent(in) :: interval_b
    real, dimension(2), intent(out) :: intervals_c
    integer, intent(out) :: n_intervals_c

    real :: min_a, max_a, min_b, max_b

    min_a = minval(interval_a)
    max_a = maxval(interval_a)
    min_b = minval(interval_b)
    max_b = maxval(interval_b)

    if(max_b <= min_a .or. min_b >= max_a) then
      n_intervals_c = 0
    else
      intervals_c(1) = max(min_a, min_b)
      intervals_c(2) = min(max_a, max_b)
      n_intervals_c = 1
    end if

  end subroutine intersect_intervals_rank_1
  
  pure subroutine intersect_intervals_rank_2(interval_a, interval_b, intervals_c, n_intervals_c)
    real, dimension(2), intent(in) :: interval_a
    real, dimension(2), intent(in) :: interval_b
    real, dimension(2, interval_buf_size), intent(out) :: intervals_c
    integer, intent(out) :: n_intervals_c
    
    call intersect_intervals(interval_a, interval_b, intervals_c(:, 1), n_intervals_c)
  
  end subroutine intersect_intervals_rank_2
  
  pure subroutine intersect_intervals_rank_3(interval_a, interval_b, intervals_c, n_intervals_c)
    real, dimension(1, 2), intent(in) :: interval_a
    real, dimension(1, 2), intent(in) :: interval_b
    real, dimension(1, 2, interval_buf_size), intent(out) :: intervals_c
    integer, intent(out) :: n_intervals_c
    
    call intersect_intervals(interval_a(1, :), interval_b(1, :), intervals_c(1, :, 1), n_intervals_c)
  
  end subroutine intersect_intervals_rank_3

end module libsupermesh_interval_intersection
