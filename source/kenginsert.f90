!This subroutine makes an array of the potential from an insertion of the Gaussian Basis
SUBROUTINE kenginsert(kinetic, gridpts, gridweights)
USE inputvariables
IMPLICIT NONE
!
INTEGER :: i, j, k, norbits
REAL(KIND=8) :: lagin(1:numgauss, 1:numgauss), gaussexp(1:numgauss)
REAL(KIND=8) :: overlaps(1:numgauss,1:numgauss), inverse_overlaps(1:numgauss,1:numgauss)
COMPLEX*16 :: lagmat(1:nbasis,1:numgauss)
COMPLEX*16 :: kinetic(1:nbasis, 1:nbasis)
COMPLEX*16 :: gridpts(1:nbasis), gridweights(1:nbasis), temp(1:nbasis,1:numgauss)
COMPLEX*16 :: smalltemp(1:numgauss,1:numgauss)
REAL(KIND=8) :: orthorbitals(1:numgauss,1:numgauss)
!
!Read in gaussian exponents, overlaps and the potential energy in Gaussians
CALL lagkenet(gaussexp, overlaps, lagin)

! CALL SVD routine to calculate the orthonormal orbitals or the pseudoinverse
CALL SVD_Ortho(overlaps, inverse_overlaps, numgauss, orthorbitals, norbits, SVD_tol)
WRITE(6,*) "SVD completed"

!Create the Gaussian matrix evaluated in dvrs with the appropriate weights and factors
CALL gaussmat(gaussexp, gridpts, gridweights, lagmat)

! For Psudeo inverse
IF(switchv .EQ. 1) THEN

   WRITE(*,*) "Using the SVD pseudoinverse for KE"
   WRITE(*,*) "SVD tolerance is: ", SVD_tol

   !Create the inserted potential evaluated on the DVR
   smalltemp = MATMUL(inverse_overlaps,MATMUL(lagin,inverse_overlaps))
   kinetic = MATMUL(lagmat,MATMUL(smalltemp,TRANSPOSE(lagmat)))

! For orthogonal orbitals
ELSE IF(switchv .EQ. 2) THEN
   WRITE(*,'(1x,"Using SVD orthonormal orbitals for KE with a tolerance of ", ES15.6)')  SVD_tol
   WRITE(*,'(1x,"There are ", I5, " orthonormal orbitals")') norbits

   ! Evaluate the insertion multiplication
   smalltemp = MATMUL(TRANSPOSE(orthorbitals),MATMUL(lagin,orthorbitals))
   smalltemp = MATMUL(orthorbitals,MATMUL(smalltemp,TRANSPOSE(orthorbitals)))
   kinetic = MATMUL(lagmat,MATMUL(smalltemp,TRANSPOSE(lagmat)))
ENDIF

ENDSUBROUTINE
