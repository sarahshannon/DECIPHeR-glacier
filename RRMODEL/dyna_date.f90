!-----------------------------------------------------------------------------
!--- Sarah Shannon 03/12/2020 date functions for the glacier model
!-----------------------------------------------------------------------------
module dyna_date

  implicit none

contains

 logical function is_leap_year(year)

   implicit none
   
!-----------------------------------------------------------------------------
! Description:
!   Returns .true. if the year is a leap year
!   Returns .false. if not
!-----------------------------------------------------------------------------

   integer, intent(in) :: year

   is_leap_year = ( mod(year, 4) == 0 .and. mod(year, 100) /= 0 ) .or.       &
        ( mod(year, 400) == 0 )

   return

 end function is_leap_year

 
 integer function days_in_year(year, l_leap)

   implicit none
   
!-----------------------------------------------------------------------------
!   Returns number of days in the year 
!-----------------------------------------------------------------------------

    integer,intent(in) :: year
  
    logical, intent(in) :: l_leap  ! true for leap years/false for no leap years


! Default number of days in year is (obviously) 365
    days_in_year = 365

    IF ( l_leap .AND. is_leap_year(year) ) days_in_year = 366

    RETURN

  END FUNCTION days_in_year


  INTEGER FUNCTION days_in_month(year, month, l_leap)

    IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Returns the number of days in the given month for the given year
!-----------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: year, month
   
    LOGICAL, INTENT(IN) :: l_leap  ! T for leap years
                                   ! F for no leap years


! If using 'normal' year then we need to do different things for each month
    SELECT CASE ( month )

      CASE ( 2 )
        IF ( l_leap .AND. is_leap_year(year) ) THEN
          days_in_month = 29
        ELSE
          days_in_month = 28
        END IF
!  months with 30 days 
      CASE ( 4, 6, 9, 11 )
        days_in_month = 30
! all others have 31 days
      CASE DEFAULT
        days_in_month = 31
    END SELECT

    RETURN

  END FUNCTION days_in_month


  INTEGER FUNCTION day_of_year(year, month, day, l_leap)

    IMPLICIT NONE

!-----------------------------------------------------------------------------
!   Returns the day of the year that corresponds to year, month and day
!-----------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: year, month, day
   
    LOGICAL, INTENT(IN) :: l_leap  ! T for leap years
                                   ! F for no leap years

! Work variables
    INTEGER :: m  ! Month counter


    day_of_year = 0
    m = 0
! Count through the full months worth of days we need to add
    DO
      m = m + 1
      IF ( m >= month ) EXIT

      day_of_year = day_of_year + days_in_month(year, m, l_leap)
    END DO

! Add the remaining days
    day_of_year = day_of_year + day

    RETURN

  END FUNCTION day_of_year

END MODULE dyna_date
