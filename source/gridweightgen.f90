!This subroutine generates the weights to normalize the basis functions.
!It's used on the right hand side of the driven equation
SUBROUTINE gridweightgen(weights, gridweights)
USE inputvariables
IMPLICIT NONE
!
INTEGER :: i, j, k
COMPLEX*16 :: weights(0:numelements-1,0:norder-1)
COMPLEX*16 :: gridweights(1:nbasis)
!
OPEN(unit=9994, file='gridweights.out', status='unknown')
!
!Generates the basis functions
!
k=1
DO j=0, numelements-1
   DO i=1, norder-2
      gridweights(j*(norder-1)+i) = weights(j,i)
   ENDDO
   IF(j .NE. numelements-1) THEN
      gridweights((j+1)*(norder-1)) = weights(j,norder-1)+weights(j+1,0)
   ENDIF
ENDDO
!
!Output the basis weights
!
DO i=1,nbasis
   WRITE(9994,*) DREAL(gridweights(i)), DIMAG(gridweights(i))
ENDDO

CLOSE(9994)
!
END SUBROUTINE
