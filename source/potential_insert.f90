! This subroutine makes an array of the potential from an insertion of the Gaussian Basis
SUBROUTINE potential_insert(potential, gridpts, gridweights, vnucf)
USE inputvariables
IMPLICIT NONE

! Input variables
COMPLEX(KIND=DBL), INTENT(IN) :: gridpts(1:nbasis)
COMPLEX(KIND=DBL), INTENT(IN) :: gridweights(1:nbasis)
INTEGER, INTENT(IN) :: vnucf ! include, exclude or only vnuc as the potential

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
COMPLEX(KIND=DBL) :: contracted(1:numgauss,1:nbasis)

! Local temps
INTEGER :: i, j, k
INTEGER :: part=0
COMPLEX(KIND=DBL) ::  temp(1:nbasis,1:numgauss)
COMPLEX(KIND=DBL) :: smalltemp(1:numgauss,1:numgauss)
COMPLEX(KIND=DBL) :: largetemp(1:numprimg,1:numprimg)

!Read in gaussian exponents, overlaps and the potential energy in Gaussians
CALL gausspot(numgauss, overlaps, gaussmatin, vnucf)

! CALL SVD routine to calculate the orthonormal orbitals or the pseudoinverse
CALL SVD_Ortho(overlaps, inverse_overlaps, numgauss, orthorbitals, norbits, SVD_tol)
WRITE(6,*) "SVD completed"

IF (vnucf .LT. 0) part = 1
!Create the Gaussian matrix evaluated in dvrs with the appropriate weights and factors
CALL gaussmat(gridpts, gridweights, gaussdvr, mfac, part)

! Formal complex conjugate of the angular functions
mfac=-1
CALL gaussmat(gridpts, gridweights, conjgaussdvr, mfac, part)

!!$! Call for contraction matrix
!!$OPEN(UNIT=95, FILE='xform.dat', STATUS='OLD', ACTION='READ')
!!$DO i=1, numgauss
!!$   READ(95, *) (xformmat(i,j),j=1,numprimg)
!!$ENDDO
!!$CLOSE(95)
!!$
!!$! Write out the contractions
!!$contracted  = MATMUL(xformmat, TRANSPOSE(gaussdvr))
!!$
!!$OPEN(UNIT=110, FILE='contraction.out', STATUS='UNKNOWN', ACTION='WRITE')
!!$DO i=1, nbasis
!!$   WRITE(110, '(60ES18.8E3)') gridpts(i), (contracted(j,i), j=1, numgauss)
!!$ENDDO
!!$CLOSE(110)

!!$OPEN(UNIT=115, FILE='gprim.out', STATUS='UNKNOWN', ACTION='WRITE')
!!$DO i=1, nbasis
!!$   WRITE(115, '(60ES18.8E3)') (gaussdvr(i,j), j=1, numgauss)
!!$ENDDO
!!$CLOSE(115)

! Use SVD orbitals or inverse for the insertion?
!    SVDswitch == 1: orbitals
!    SVDswitch == 0: inverse
insertinverse: IF(switchv .EQ. 2) THEN
   WRITE(*,'(1x,"Using SVD orthonormal orbitals for the potential with a tolerance of ", ES15.6)')  SVD_tol
   WRITE(*,'(1x,"There are ", I5, " orthonormal orbitals")') norbits

   ! Evaluate the insertion multiplication
   smalltemp = MATMUL(TRANSPOSE(orthorbitals),MATMUL(gaussmatin,orthorbitals))
   smalltemp = MATMUL(orthorbitals,MATMUL(smalltemp,TRANSPOSE(orthorbitals)))
   potential = MATMUL(conjgaussdvr,MATMUL(smalltemp,TRANSPOSE(gaussdvr)))


ELSE IF(switchv .EQ. 1) THEN insertinverse
   WRITE(*,*) "Using the SVD pseudoinverse for the potential"
   WRITE(*,*) "SVD tolerance is: ", SVD_tol

   !Create the inserted potential evaluated on the DVR
   smalltemp = MATMUL(inverse_overlaps,MATMUL(gaussmatin,inverse_overlaps))
   potential = MATMUL(conjgaussdvr,MATMUL(smalltemp,TRANSPOSE(gaussdvr)))


ELSE insertinverse
   WRITE(*,*) "Insertion potential is undefined"
   STOP
ENDIF insertinverse

! Insertion sucessful!
WRITE(*,*) "Potential Energy Inserted"

ENDSUBROUTINE