!This subroutine evaluates and diagonalizes the projector from a Gaussian basis
SUBROUTINE projection(gridpts, gridweights)
USE inputvariables
IMPLICIT NONE
!
INTEGER :: i, j, k, info
CHARACTER(LEN=1) :: jobvl = 'N'
REAL(KIND=8) :: gaussexp(1:numgauss)
REAL(KIND=8) :: overlaps(1:numgauss,1:numgauss)
COMPLEX*16 :: lagmat(1:nbasis,1:numgauss), lagmattran(1:numgauss,1:nbasis)
COMPLEX*16 :: projector(1:nbasis, 1:nbasis)
COMPLEX*16 :: gridpts(1:nbasis), gridweights(1:nbasis), temp(1:nbasis,1:numgauss)
COMPLEX*16 :: eigen(1:nbasis), empt(1:nbasis,1:nbasis)
COMPLEX*16 :: workspace(1:20*nbasis)
REAL(KIND=8) :: cwork(1:2*nbasis)
!
!Read in the Gaussian exponents
OPEN(unit=1, file='gexp.dat', status='UNKNOWN', action='READ')
DO i=1, numgauss
    READ(1,'(1x,80ES17.9)') gaussexp(i) 
ENDDO
CLOSE(1)
!
!Read in the overlaps of the Gaussian exponents
OPEN(unit=2, file='overlaps.dat', status='UNKNOWN',action='READ')
DO i=1, numgauss
    READ(2,'(1x,80ES17.9)') (overlaps(i,j), j=1, numgauss)  
ENDDO
CLOSE(2)
!
!Invert the overlap matrix
CALL SmxInv(overlaps, numgauss)
!
!Evaluate the Guassians on the DVR with the correct factors
CALL gaussmat(gaussexp, gridpts, gridweights, lagmat)
!
!Create the projector matrix on the DVR
temp = MATMUL(lagmat,overlaps)
lagmattran = TRANSPOSE(lagmat)
projector = MATMUL(temp,lagmattran)
WRITE(*,*) "Projection done!"
!
!Diagonalize the projector
CALL ZGEEV(jobvl, jobvl, nbasis, projector, nbasis, eigen, empt, nbasis, empt, nbasis, workspace, 20*nbasis, cwork, info)
!
WRITE(*,*) "Did zgeev work?", info
!
OPEN(unit=5555, file='proj.dat', status='UNKNOWN',action='WRITE')
DO i=1, nbasis
   WRITE(5555,'(1x,80ES17.9)') DREAL(eigen(i)), DIMAG(eigen(i))
ENDDO
CLOSE(5555)
!
ENDSUBROUTINE
