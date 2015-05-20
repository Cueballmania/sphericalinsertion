SUBROUTINE coef(sym, expo, norm)
! This subroutine returns the normalization coeffecient for Mesa
! Cartesian gaussians.
IMPLICIT NONE

! Declare double precision
INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=13,r=200)
REAL(KIND=DBL), PARAMETER :: PI = 3.141592653589

! Input variables
CHARACTER(LEN=2), INTENT(IN) :: sym
REAL(KIND=DBL), INTENT(IN) :: expo

! Output variable
REAL(KIND=DBL), INTENT(OUT) :: norm

IF(sym .EQ. 'ss') THEN
   norm = (2.0d0*expo/PI)**(0.75) 
ELSEIF(sym .EQ. 'px' .OR. sym .EQ. 'py' .OR. sym .EQ. 'pz') THEN
   norm = (2.0d0**7*expo**5/PI**3)**(0.25)
ELSEIF(sym .EQ. 'xx' .OR. sym .EQ. 'yy' .OR. sym .EQ. 'zz') THEN
   norm = (2.0d0**11*expo**7/PI**3)**(0.25)
ELSEIF(sym .EQ. 'xy' .OR. sym .EQ. 'xz' .OR. sym .EQ. 'yz') THEN
   norm = (2.0d0**11*expo**7/PI**3)**(0.25)
ELSE
   WRITE(*,*) "Sym = ", sym, " not recognized for expo=", expo
ENDIF

END SUBROUTINE
