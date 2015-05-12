!This subroutine makes an array of the potential from Laguerre basis
!The potentail dies off to zero on the complex scale
SUBROUTINE coulpot_sch2(lagpotential, gridpts, gridweights)
USE inputvariables
IMPLICIT NONE
!
INTEGER :: i, j,k
REAL(KIND=8) :: lagin(1:numgauss, 1:numgauss), gaussexp(1:numgauss), overlaps(1:numgauss, 1:numgauss)
COMPLEX*16 :: lagmat(1:nbasis,1:numgauss), lagmattran(1:numgauss,1:nbasis)
COMPLEX*16 :: lagpotential(1:nbasis, 1:nbasis)
COMPLEX*16 :: gridpts(1:nbasis), gridweights(1:nbasis), temp(1:nbasis,1:numgauss)
!
CALL lagpot(gaussexp, overlaps, lagin)
!
CALL SmxInv(lagin, numgauss)
!
OPEN(unit=1005, file='lagpotent.out', status='unknown')
!
CALL gaussmat_sch(gaussexp, overlaps, gridpts, lagmat)
!
temp = MATMUL(lagmat,lagin)
!
lagmattran = TRANSPOSE(lagmat)
!
WRITE(*,*) "TRANS done"
lagpotential = MATMUL(temp,lagmattran)
WRITE(*,*) "V done!"
!
DO i=1, nbasis
   DO j=1, nbasis
      lagpotential(i,j) = lagpotential(i,j)*ZSQRT(gridweights(i)*gridweights(j))
   ENDDO
   WRITE(1005,*) (lagpotential(i,k),k=1,nbasis)
ENDDO
CLOSE(1005)
!
ENDSUBROUTINE
