! This function calculates the value of the partitioning function at a given x between [a,b]
FUNCTION partitioning(x, a, b)
IMPLICIT NONE

! Declare double precision
INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=13,r=200)
REAL(KIND=DBL) :: partitioning

! Input variables: value and the boundaries a<b
REAL(KIND=DBL), INTENT(IN) :: x
REAL(KIND=DBL), INTENT(IN) :: a, b

REAL(KIND=DBL) :: xtrans

xtrans = (x-a)/(b-a) - 1.0d0

partitioning = EXP(xtrans*xtrans/(xtrans*xtrans-1.0d0))
RETURN
END FUNCTION partitioning
