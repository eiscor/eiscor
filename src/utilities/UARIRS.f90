#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! UARIRS (Utility Auxiliary Routine Initialize Random Seed)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine initializes the random number generator using the cpu
! clock.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
! OUTPUT VARIABLES:
!
! INFO            INTEGER
!                    INFO = 1 implies array allocation failed
!                    INFO = 0 implies successful computation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine UARIRS(INFO)

  implicit none
  
  ! input variables
  integer, intent(inout) :: INFO

  ! compute variables
  integer :: ii, n, clock
  integer, allocatable :: seed(:)
  
  ! initialize INFO
  INFO = 0
  
  ! get size of see        
  call random_seed(size = n)
  
  ! allocate memory for seed
  allocate(seed(n))
  
  ! check allocation
  if (allocated(seed).EQV..FALSE.) then
    INFO = 1
  
    ! print error in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"Array allocation failed",INFO,INFO)
    end if 
    
    return
  end if 
  
  ! get current clock time        
  call system_clock(count=clock)
  
  ! store seeds        
  seed = clock + 37 * (/ (ii - 1, ii = 1, n) /)
  
  ! set the generator
  call random_seed(put = seed)
  
  ! free memory        
  deallocate(seed)

end subroutine UARIRS
