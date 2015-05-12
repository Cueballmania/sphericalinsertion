SUBROUTINE quadprime(felement, braval, ketval, integ, weights)
USE inputvariables
IMPLICIT NONE
!
INTEGER :: i, j, felement
COMPLEX*16 :: integ
COMPLEX*16 :: braval(0:norder-1), ketval(0:norder-1)
COMPLEX*16 :: weights(0:numelements-1,0:norder-1)
!
integ = czero
DO i=0, norder-1
   integ = integ + braval(i)*ketval(i)*weights(felement,i)
ENDDO
!
IF(felement.GT.relements) integ=integ/contourphase
!
ENDSUBROUTINE
