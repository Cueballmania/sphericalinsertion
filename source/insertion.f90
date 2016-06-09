! This subroutine makes an array of the potential from an insertion of the Gaussian Basis
SUBROUTINE insertion(potential, gridpts, gridweights)
USE inputvariables
IMPLICIT NONE

! Input variables
COMPLEX(KIND=DBL), INTENT(IN) :: gridpts(1:nbasis)
COMPLEX(KIND=DBL), INTENT(IN) :: gridweights(1:nbasis)

! Output potential
COMPLEX(KIND=DBL), INTENT(OUT) :: potential(1:nbasis, 1:nbasis)

! Read-in matrices from mesaread
REAL(KIND=DBL) :: overlaps(1:numgauss,1:numgauss)
REAL(KIND=DBL) :: gaussmatin(1:numgauss, 1:numgauss)

! Orthogonalization matrices
INTEGER :: norbits
REAL(KIND=DBL) :: inverse_overlaps(1:numgauss,1:numgauss)
REAL(KIND=DBL) :: orthorbitals(1:numgauss,1:numgauss)

! DVR/Gaussian xform
INTEGER :: mfac = 1
COMPLEX(KIND=DBL) :: gaussdvr(1:nbasis,1:numprimg)
COMPLEX(KIND=DBL) :: conjgaussdvr(1:nbasis,1:numprimg)

REAL(KIND=DBL) :: xformmat(1:numgauss,1:numprimg)

! Local temps
INTEGER :: i, j, k
INTEGER :: part=0
COMPLEX(KIND=DBL) ::  temp(1:nbasis,1:numgauss)
COMPLEX(KIND=DBL) :: smalltemp(1:numgauss,1:numgauss)
COMPLEX(KIND=DBL) :: largetemp(1:numprimg,1:numprimg)

!Read in gaussian exponents, overlaps and the potential energy in Gaussians
CALL gausspot(numgauss, overlaps, gaussmatin, kedvr, switchv, partitionflag)

! CALL SVD routine to calculate the orthonormal orbitals or the pseudoinverse
CALL SVD_Ortho(overlaps, inverse_overlaps, numgauss, orthorbitals, norbits, SVD_tol)
WRITE(6,*) "SVD completed"

!Create the Gaussian matrix evaluated in dvrs with the appropriate weights and factors
CALL gaussmat(gridpts, gridweights, gaussdvr, mfac)

! Formal complex conjugate of the angular functions
mfac=-1
CALL gaussmat(gridpts, gridweights, conjgaussdvr, mfac)

IF (numgauss .EQ. numprimg) THEN
   xformmat = 0.0d0
   DO i=1, numgauss
      xformmat(i,i) = 1.0d0
   ENDDO
   
ELSE
   ! Call for contraction matrix
   OPEN(UNIT=95, FILE='xform.dat', STATUS='OLD', ACTION='READ')
   DO i=1, numgauss
      READ(95, *) (xformmat(i,j),j=1,numprimg)
   ENDDO
   CLOSE(95)
ENDIF


! Use SVD orbitals or inverse for the insertion?
!    SVDswitch == 1: orbitals
!    SVDswitch == 0: inverse
insertinverse: IF(switchv .GE. 2) THEN
   WRITE(*,'(1x,"Using SVD orthonormal orbitals for the potential with a tolerance of ", ES15.6)')  SVD_tol
   WRITE(*,'(1x,"There are ", I5, " orthonormal orbitals")') norbits

   ! Evaluate the insertion multiplication
   smalltemp = MATMUL(TRANSPOSE(orthorbitals),MATMUL(gaussmatin,orthorbitals))
   smalltemp = MATMUL(orthorbitals,MATMUL(smalltemp,TRANSPOSE(orthorbitals)))
   largetemp = MATMUL(TRANSPOSE(xformmat),MATMUL(smalltemp,xformmat))
   potential = MATMUL(conjgaussdvr,MATMUL(largetemp,TRANSPOSE(gaussdvr)))


ELSE IF(switchv .EQ. 1) THEN insertinverse
   WRITE(*,*) "Using the SVD pseudoinverse for the potential"
   WRITE(*,*) "SVD tolerance is: ", SVD_tol

   !Create the inserted potential evaluated on the DVR
   smalltemp = MATMUL(inverse_overlaps,MATMUL(gaussmatin,inverse_overlaps))
   largetemp = MATMUL(TRANSPOSE(xformmat),MATMUL(smalltemp,xformmat))
   potential = MATMUL(conjgaussdvr,MATMUL(largetemp,TRANSPOSE(gaussdvr)))


ELSE insertinverse
   WRITE(*,*) "Insertion potential is undefined"
   STOP
ENDIF insertinverse

! Insertion sucessful!
WRITE(*,*) "Potential Energy Inserted"

ENDSUBROUTINE
