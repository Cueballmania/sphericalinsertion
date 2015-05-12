SUBROUTINE calcheck(info)
! This subroutine checks that the calculation is self consistant
! and all of the input parameters make sense
! This subroutine also output calculation information to the screen.
! Returns an integer value of zero if the input is okay
USE inputvariables
IMPLICIT NONE

! Output variable
INTEGER, INTENT(OUT) :: info

info=0

! Check finite elements
IF(celements > numelements) info = 1 ! There are too many complex elements
IF(angle > 45.0d0) info = 1 ! Angle for scaling is greater than Pi/4 = 45 degress
IF(angle < 0.0) info = 1 ! Negative scaling angle

WRITE(*,'(" There are ", I3, " elements.  The last ", I3," are complex")') numelements, celements
IF(celements > 0) THEN
   WRITE(*, '("The angle of scaling is ", F8.3, " degrees or ", F8.6, " pi radians")') angle, angle*Pi/180.0
ENDIF

! Check if the number of quadrature points per eleement less than MAXN
IF(norder > MAXN) info = 2 ! There are too many points!

WRITE(*,'(" There are ", I4," quadrature points per element giving ", I5, " total points")') norder, nbasis

! Check the angular quantum numbers
IF(l_ang < 0) info = 3 ! L cannot be a negative number
IF(ABS(m_ang) > l_ang) info = 3 ! M > L which is forbidden

WRITE(*,'(" The calculation uses L=", I2," and M =",I2)') l_ang, m_ang

! Specify the type of calculation
WRITE(*,*) "The calculation types are"
IF(kedvr .EQ. 0)  THEN
   WRITE(*,*) "DVR Kinetic Energy"
ELSE
   WRITE(*,*) "Insertion Kinetic Energy"
ENDIF

IF(switchv .EQ. 0) THEN
   WRITE(*,*) "switchv=0: DVR Potential Energy"
ELSEIF(switchv .EQ. 1) THEN
   WRITE(*,*) "switchv=1: Gaussian Potential Insertion (Psudeoinverse)"
ELSEIF(switchv .EQ. 2) THEN
   WRITE(*,*) "switchv=2: Gaussian Potential Insertion (Orthogonal Orbitals)"
ENDIF


! Check the gaussian insertion
IF(numgauss > numprimg .AND. (switchv /= 0 .OR. kedvr /= 0)) info = 4 ! More contracted Guassians than primative

IF(switchv /= 0 .OR. kedvr /= 0) THEN
   WRITE(*,'("The insertion calculation has ", I3, " contracted Guassians and ", I3, " primatives")') numgauss, numprimg
ENDIF

END SUBROUTINE
