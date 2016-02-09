#include "fdebug.h"

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

  pure subroutine intersect_intervals_rank_1(intervalA, intervalB, intervalC, n_intervalsC)
    real, dimension(2), intent(in) :: intervalA
    real, dimension(2), intent(in) :: intervalB
    real, dimension(2), intent(out) :: intervalC
    integer, intent(out) :: n_intervalsC

    real :: min_a, max_a, min_b, max_b

    min_a = minval(intervalA)
    max_a = maxval(intervalA)
    min_b = minval(intervalB)
    max_b = maxval(intervalB)

    if(max_b <= min_a .or. min_b >= max_a) then
      n_intervalsC = 0
    else
      intervalC(1) = max(min_a, min_b)
      intervalC(2) = min(max_a, max_b)
      n_intervalsC = 1
    end if

  end subroutine intersect_intervals_rank_1
  
  pure subroutine intersect_intervals_rank_2(intervalA, intervalB, intervalsC, n_intervalsC)
    real, dimension(2), intent(in) :: intervalA
    real, dimension(2), intent(in) :: intervalB
    real, dimension(2, interval_buf_size), intent(out) :: intervalsC
    integer, intent(out) :: n_intervalsC
    
    call intersect_intervals(intervalA, intervalB, intervalsC(:, 1), n_intervalsC)
  
  end subroutine intersect_intervals_rank_2
  
  pure subroutine intersect_intervals_rank_3(intervalA, intervalB, intervalsC, n_intervalsC)
    real, dimension(1, 2), intent(in) :: intervalA
    real, dimension(1, 2), intent(in) :: intervalB
    real, dimension(1, 2, interval_buf_size), intent(out) :: intervalsC
    integer, intent(out) :: n_intervalsC
    
    call intersect_intervals(intervalA(1, :), intervalB(1, :), intervalsC(1, :, 1), n_intervalsC)
  
  end subroutine intersect_intervals_rank_3

end module libsupermesh_interval_intersection
