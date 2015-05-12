SUBROUTINE eigentest(kinetic,potential, grid, weights, wavefunction)
! This subroutine tests various things about the eigenvariables and the Hamiltonian
USE inputvariables
IMPLICIT NONE

! Input variables
COMPLEX(KIND=DBL), INTENT(IN) :: kinetic(1:nbasis,1:nbasis)
COMPLEX(KIND=DBL), INTENT(IN) :: potential(1:nbasis,1:nbasis)
COMPLEX(KIND=DBL), INTENT(IN) :: wavefunction(1:nbasis)
COMPLEX(KIND=DBL), INTENT(IN) :: grid(1:nbasis)
COMPLEX(KIND=DBL), INTENT(IN) :: weights(1:nbasis)

! Local matrices
COMPLEX(KIND=DBL) :: ham(1:nbasis,1:nbasis)
COMPLEX(KIND=DBL) :: gaussdvr(1:nbasis,1:numprimg)
COMPLEX(KIND=DBL) :: conjgaussdvr(1:nbasis,1:numprimg)
REAL(KIND=DBL) :: xformmat(1:numgauss,1:numprimg)
REAL(KIND=DBL) :: mesaeigen(1:numgauss,1:numgauss)
COMPLEX(KIND=DBL) :: eigens(1:numgauss,1:numgauss)
REAL(KIND=DBL) :: overlaps(1:numgauss,1:numgauss)
COMPLEX(KIND=DBL) :: dvrovers(1:nbasis,1:nbasis)
COMPLEX(KIND=DBL) :: scaled_ham(1:numgauss,1:numgauss)
COMPLEX(KIND=DBL) :: newoverlap(1:numgauss,1:numgauss)
COMPLEX(KIND=DBL) :: gstate(1:nbasis)

! Local variables
INTEGER :: i,j
INTEGER :: mfac = 1

! LAPACK variables
INTEGER :: info
INTEGER :: lwork
REAL(KIND=DBL) :: rwork(1:8*numgauss)
COMPLEX(KIND=DBL) :: work(1:2*numgauss)
COMPLEX(KIND=DBL), DIMENSION(1:nbasis) :: alpha, beta
COMPLEX(KIND=DBL), DIMENSION(1:nbasis,1:nbasis) :: vl, vr
CHARACTER(LEN=1) :: jobvl='N', jobvr='N'

lwork=numgauss*2

ham=kinetic+potential

! Read in the eigenfunctions from the Mesa Calculation
OPEN(UNIT=15, FILE='mesafunct.out', STATUS='OLD', ACTION='READ')
DO i=1, numgauss
   READ(15,*) (mesaeigen(i,j),j=1,numgauss)
ENDDO
CLOSE(15)

! Call for contraction matrix
OPEN(UNIT=95, FILE='xform.dat', STATUS='OLD', ACTION='READ')
DO i=1, numgauss
   READ(95, *) (xformmat(i,j),j=1,numprimg)
ENDDO
CLOSE(95)

! Generate the dvr/gaussian overlaps
CALL gaussmat(grid, weights, gaussdvr, mfac)
mfac=-1
CALL gaussmat(grid, weights, conjgaussdvr, mfac)

! Create the ground state Mesa eigenfunction
gstate = MATMUL(mesaeigen(:,1), MATMUL(xformmat, TRANSPOSE(conjgaussdvr)))

! Transform the hamiltonian 
scaled_ham = MATMUL(xformmat,MATMUL(TRANSPOSE(conjgaussdvr),MATMUL(ham,MATMUL(gaussdvr,TRANSPOSE(xformmat)))))

eigens = MATMUL(TRANSPOSE(mesaeigen),MATMUL(scaled_ham,mesaeigen))

OPEN(UNIT=1212,FILE='eigen.out',STATUS='UNKNOWN',ACTION='WRITE')
DO i=1, numgauss
   WRITE(1212,'(80ES16.6)') (eigens(i,j), j=1, numgauss)
ENDDO
CLOSE(1212)

OPEN(UNIT=1000,FILE='eigendiag.out',STATUS='UNKNOWN',ACTION='WRITE')
DO i=1, numgauss
   WRITE(1000,'(80ES16.6)') eigens(i,i)
ENDDO
CLOSE(1000)

! Read in overlaps
CALL readmesa('overlaps.dat', numgauss, overlaps)

OPEN(UNIT=9999,FILE='wavefunction.out',STATUS='UNKNOWN',ACTION='WRITE')
DO i=1, nbasis
   WRITE(9999,'(6ES18.8E3)') REAL(grid(i)), AIMAG(grid(i)) &
        ,REAL(wavefunction(i)/SQRT(weights(i))), AIMAG(wavefunction(i)/SQRT(weights(i)))&
        ,REAL(gstate(i)/SQRT(weights(i))), AIMAG(gstate(i)/SQRT(weights(i)))
ENDDO
CLOSE(9999)

dvrovers = MATMUL(conjgaussdvr,MATMUL(TRANSPOSE(xformmat),MATMUL(overlaps,MATMUL(xformmat,TRANSPOSE(gaussdvr)))))

!!$CALL ZGGEV(jobvl, jobvr, numgauss, ham, numgauss, dvrovers, numgauss, alpha, beta, vl, numgauss, vr, &
!!$     & numgauss, work, lwork, rwork, info)
!!$
!!$DO i=1, numgauss
!!$   WRITE(120,*) alpha(i), beta(i)
!!$ENDDO

! Check the overlap matrix - compute the difference and find the eigenvalues
newoverlap = MATMUL(xformmat,MATMUL(MATMUL(TRANSPOSE(conjgaussdvr),gaussdvr),TRANSPOSE(xformmat)))

OPEN(UNIT=104,FILE='overdiff.out',STATUS='UNKNOWN',ACTION='WRITE')
DO i=1, numgauss
   WRITE(104,'(100ES16.6)') (DREAL(newoverlap(i,j)),j=1,numgauss)
ENDDO
CLOSE(104)

ENDSUBROUTINE
