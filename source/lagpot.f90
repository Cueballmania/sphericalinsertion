!This subroutine constructs the matrix of the Laguerre potential energy
!from the old code's lower diagonal matrix
SUBROUTINE lagpot(gaussexp, overlaps, matrixin)
USE inputvariables
IMPLICIT NONE
!
INTEGER :: i, j
REAL*8 :: matrixin(1:numgauss, 1:numgauss), gaussexp(1:numgauss)
REAL*8 :: overlaps(1:numgauss, 1:numgauss)
!
!Read in exponents
OPEN(unit=1, file='gexp.dat', status='UNKNOWN', action='READ')
DO i=1, numgauss
    READ(1,'(1x,80ES17.9)') gaussexp(i) 
ENDDO
CLOSE(1)
!
!Read in overlaps
OPEN(unit=2, file='overlaps.dat', status='UNKNOWN',action='READ')
DO i=1, numgauss
    READ(2,'(1x,80ES17.9)') (overlaps(i,j), j=1, numgauss)  
ENDDO
CLOSE(2)
!
!Read in potential
OPEN(unit=3, file='gpot.dat', status='UNKNOWN',action='READ')
DO i=1, numgauss
    READ(3,'(1x,80ES17.9)') (matrixin(i,j), j=1, numgauss)  
ENDDO
CLOSE(3)
!
ENDSUBROUTINE
