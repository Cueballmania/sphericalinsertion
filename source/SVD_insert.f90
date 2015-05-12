!This subroutine makes an array of the potential from an insertion of the Gaussian Basis
SUBROUTINE svd_insert(lagpotential, gridpts, gridweights)
USE inputvariables
IMPLICIT NONE
!
INTEGER :: i, j, k, norbits
INTEGER :: mfac = 1
REAL(KIND=DBL) :: lagin(1:numgauss, 1:numgauss), gaussexp(1:numgauss)
REAL(KIND=DBL) :: overlaps(1:numgauss,1:numgauss), inverse_overlaps(1:numgauss,1:numgauss)
COMPLEX(KIND=DBL) :: gaussdvr(1:nbasis,1:numprimg)
COMPLEX(KIND=DBL) :: conjgaussdvr(1:nbasis,1:numprimg)
COMPLEX(KIND=DBL) :: lagpotential(1:nbasis, 1:nbasis)
COMPLEX(KIND=DBL) :: gridpts(1:nbasis), gridweights(1:nbasis), temp(1:nbasis,1:numgauss)
REAL(KIND=DBL) :: xformmat(1:numgauss,1:numprimg)
REAL(KIND=DBL) :: orthorbitals(1:numgauss,1:numgauss)
COMPLEX(KIND=DBL) :: smalltemp(1:numgauss,1:numgauss)
COMPLEX(KIND=DBL) :: largetemp(1:numprimg,1:numprimg)

COMPLEX(KIND=DBL) :: newoverlap(1:numgauss,1:numgauss)
REAL(KIND=DBL):: overdiff(1:numgauss,1:numgauss)

!Read in gaussian exponents, overlaps and the potential energy in Gaussians
CALL svd_lagpot(overlaps, lagin)

! CALL SVD routine to calculate the orthonormal orbitals or the pseudoinverse
CALL SVD_Ortho(overlaps, inverse_overlaps, numgauss, orthorbitals, norbits, SVD_tol)
WRITE(6,*) "SVD completed"

!Create the Gaussian matrix evaluated in dvrs with the appropriate weights and factors
CALL gaussmat(gridpts, gridweights, gaussdvr, mfac)

! Formal complex conjugate of the angular functions
mfac=-1
CALL gaussmat(gridpts, gridweights, conjgaussdvr, mfac)

! Call for contraction matrix
OPEN(UNIT=95, FILE='xform.dat', STATUS='OLD', ACTION='READ')
DO i=1, numgauss
   READ(95, *) (xformmat(i,j),j=1,numprimg)
ENDDO
CLOSE(95)

! Check the overlap matrix - compute the difference and find the eigenvalues
newoverlap = MATMUL(xformmat,MATMUL(MATMUL(TRANSPOSE(conjgaussdvr),gaussdvr),TRANSPOSE(xformmat)))

OPEN(UNIT=104,FILE='overdiff.out',STATUS='UNKNOWN',ACTION='WRITE')
DO i=1, numgauss
   WRITE(104,'(100ES16.6)') (DREAL(newoverlap(i,j)),j=1,numgauss)
ENDDO
CLOSE(104)

! Use SVD orbitals or inverse for the insertion?
!    SVDswitch == 1: orbitals
!    SVDswitch == 0: inverse
insertinverse: IF(switchv .EQ. 2) THEN
   WRITE(*,'(1x,"Using SVD orthonormal orbitals for the potential with a tolerance of ", ES15.6)')  SVD_tol
   WRITE(*,'(1x,"There are ", I5, " orthonormal orbitals")') norbits

   ! Evaluate the insertion multiplication
   smalltemp = MATMUL(TRANSPOSE(orthorbitals),MATMUL(lagin,orthorbitals))
   smalltemp = MATMUL(orthorbitals,MATMUL(smalltemp,TRANSPOSE(orthorbitals)))
   largetemp = MATMUL(TRANSPOSE(xformmat),MATMUL(smalltemp,xformmat))
   lagpotential = MATMUL(conjgaussdvr,MATMUL(largetemp,TRANSPOSE(gaussdvr)))

!!$   ! Write out the orthonormal orbitals
!!$   OPEN(UNIT=1015, FILE='orthorbitals.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)  
!!$   IF(ierror == 0) THEN
!!$      DO i=1, numgauss
!!$         WRITE(1015,'(1x,800ES17.9)') (orthorbitals(i,k),k=1,numgauss)
!!$      ENDDO
!!$   ELSE
!!$      WRITE(6,*) " Error opening orthorbital file to write, ierror = ", ierror
!!$   ENDIF
!!$   CLOSE(1015)

ELSE IF(switchv .EQ. 1) THEN insertinverse
   WRITE(*,*) "Using the SVD pseudoinverse for the potential"
   WRITE(*,*) "SVD tolerance is: ", SVD_tol

   !Create the inserted potential evaluated on the DVR
   smalltemp = MATMUL(inverse_overlaps,MATMUL(lagin,inverse_overlaps))
   largetemp = MATMUL(TRANSPOSE(xformmat),MATMUL(smalltemp,xformmat))
   lagpotential = MATMUL(conjgaussdvr,MATMUL(largetemp,TRANSPOSE(gaussdvr)))


ELSE insertinverse
   WRITE(*,*) "Insertion potential is undefined"
   STOP
ENDIF insertinverse

! Insertion sucessful!
WRITE(*,*) "Potential Energy Inserted"

ENDSUBROUTINE
