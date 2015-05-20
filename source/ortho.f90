SUBROUTINE ortho(sym1, sym2, oswitch, coef, rfac)
! This subroutine determins orthogonality between two functions on the same center
! It returns the factor of their overlap divided by pi and the number of r's.
IMPLICIT NONE
  
! Declare double precision
INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=13,r=200)

! Input variables
CHARACTER(LEN=2), INTENT(IN) :: sym1, sym2

! Output variables
INTEGER, INTENT(OUT) :: oswitch
REAL(KIND=DBL), INTENT(OUT) :: coef
INTEGER, INTENT(OUT) :: rfac

! If the two are the same symmetry
IF(sym1 .EQ. sym2) THEN
   oswitch = 1
   IF(sym1 .EQ. 'ss') THEN
      coef = 4.0d0
      rfac = 0
   ELSEIF(sym1 .EQ. 'px' .OR. sym1 .EQ. 'py' .OR. sym1 .EQ. 'pz') THEN
      coef = (4.0d0/3.0d0)
      rfac = 2
   ELSEIF(sym1 .EQ. 'xx' .OR. sym1 .EQ. 'yy' .OR. sym1 .EQ. 'zz') THEN
      coef = 0.8d0
      rfac = 4
   ELSEIF(sym1 .EQ. 'xy' .OR. sym1 .EQ. 'yz' .OR. sym1 .EQ. 'xz') THEN
      coef = (4.0d0/15.0d0)
      rfac = 4
   ENDIF
ELSE
   IF(sym1 .EQ. 'ss' .AND. (sym2 .EQ. 'xx' .OR. sym2 .EQ. 'yy' .OR. sym2 .EQ. 'zz')) THEN
      oswitch = 1
      coef = (4.0d0/3.0d0)
      rfac = 2
   ELSEIF(sym2 .EQ. 'ss' .AND. (sym1 .EQ. 'xx' .OR. sym1 .EQ. 'yy' .OR. sym1 .EQ. 'zz')) THEN
      oswitch = 1
      coef = (4.0d0/3.0d0)
      rfac = 2
   ELSE
      oswitch = 0
      coef = 0.0
      rfac = 0
   ENDIF
ENDIF

END SUBROUTINE
