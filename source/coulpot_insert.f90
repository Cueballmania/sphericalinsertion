!This subroutine makes an array of the potential from an insertion of the Gaussian Basis
SUBROUTINE coulpot_insert(lagpotential, gridpts, gridweights)
USE inputvariables
IMPLICIT NONE
!
INTEGER :: i, j,k
REAL(KIND=8) :: lagin(1:numgauss, 1:numgauss), gaussexp(1:numgauss)
REAL(KIND=8) :: overlaps(1:numgauss,1:numgauss), overlapin(1:numgauss,1:numgauss)
REAL(KIND=8) :: isunit(1:numgauss,1:numgauss)
COMPLEX*16 :: lagmat(1:nbasis,1:numgauss), lagmattran(1:numgauss,1:nbasis)
COMPLEX*16 :: lagpotential(1:nbasis, 1:nbasis)
COMPLEX*16 :: gridpts(1:nbasis), gridweights(1:nbasis), temp(1:nbasis,1:numgauss)
!
!Read in gaussian exponents, overlaps and the potential energy in Gaussians
CALL lagpot(gaussexp, overlaps, lagin)
!
overlapin = overlaps
!
!Invert the overlap matrix to make the correct projector
CALL SmxInv(overlapin, numgauss)
!
isunit= MATMUL(overlaps,overlapin)
!
OPEN(unit=1015, file='isunit.out', status='unknown')
DO i=1, nbasis
   WRITE(1015,'(1x,80ES17.9)') (isunit(i,k),k=1,nbasis)
ENDDO
CLOSE(1015)
!
!Create the Gaussian matrix evaluated in dvrs with the appropriate weights and factors
CALL gaussmat(gaussexp, gridpts, gridweights, lagmat)
!
!Evaluate the insertion multiplication
temp = MATMUL(lagmat,overlaps)
temp = MATMUL(temp,lagin)
temp = MATMUL(temp,overlaps)
lagmattran = TRANSPOSE(lagmat)
lagpotential = MATMUL(temp,lagmattran)
!
WRITE(*,*) "Potential Energy Inserted"
!
OPEN(unit=1005, file='lagpotent.out', status='unknown')
DO i=1, nbasis
   WRITE(1005,'(1x,80ES17.9)') (lagpotential(i,k),k=1,nbasis)
ENDDO
CLOSE(1005)
!
ENDSUBROUTINE
