MODULE inputvariables
IMPLICIT NONE
SAVE

! Declare double precision
INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=13,r=200)

! Maximum number of quadrature points per element
INTEGER, PARAMETER :: MAXN = 500

! Static parameters
REAL(KIND=DBL), PARAMETER :: PI = 3.141592653589
REAL(KIND=DBL), PARAMETER ::TOLERANCE = 1D-10
REAL(KIND=DBL), PARAMETER :: duey=2.0d0
COMPLEX(KIND=DBL), PARAMETER :: czero=(0.0D0,0.0D0)
COMPLEX(KIND=DBL), PARAMETER :: cone=(1.0D0,0.0D0)
COMPLEX(KIND=DBL), PARAMETER :: cimag=(0.0D0,1.0D0)

! Quadrature gobals
INTEGER :: numelements, celements, relements
INTEGER :: norder, angular

! DVR gobals
INTEGER :: nbasis
REAL(KIND=DBL) :: gridend
COMPLEX(KIND=DBL) :: contourphase

! Insertion gobals
INTEGER :: numgauss, numprimg

! Calculation choices
INTEGER ::switchv, kedvr
INTEGER :: l_ang, m_ang
INTEGER :: partitionflag
REAL(KIND=DBL) :: SVD_tol
REAL(KIND=DBL) :: contourpt
REAL(KIND=DBL) :: angle
REAL(KIND=DBL) :: mass=1.0
REAL(KIND=DBL) :: outstep
END MODULE



PROGRAM mainscatter
USE inputvariables
IMPLICIT NONE

! Local variables
INTEGER :: i, j, info

! Hamiltonian matrices
COMPLEX(KIND=DBL), ALLOCATABLE :: kinetic(:,:), potential(:,:), temppot(:,:)
COMPLEX(KIND=DBL), ALLOCATABLE :: hamiltonian(:,:)

! Basis information
COMPLEX(KIND=DBL), ALLOCATABLE :: gridpts(:)
COMPLEX(KIND=DBL), ALLOCATABLE :: gridweights(:)

! Wavefunction test of properites of the ground state
INTEGER :: wfn
COMPLEX(KIND=DBL), ALLOCATABLE :: wavefunction(:)

! Partition function
REAL(KIND=DBL) :: partitioning


!The order of operations
!Read in from file
OPEN(UNIT=1000, FILE='scattering.inp', STATUS='old', ACTION='READ')

READ(1000,*) numelements, norder
READ(1000,*) l_ang, m_ang
READ(1000,*) angle, celements
READ(1000,*) numgauss, numprimg
READ(1000,*) kedvr, switchv
READ(1000,*) partitionflag
READ(1000,*) SVD_tol
READ(1000,*) wfn

CLOSE(1000)

! Calculate some constants
relements = numelements-celements
nbasis = (norder-1)*numelements-1
contourphase = DCMPLX(COS(angle*PI/180.),SIN(angle*PI/180.))

! Check input parameters and echo input parameters
CALL calcheck(info)
IF(info /= 0) THEN
   WRITE(*,'(" Calculation check failed.  info = ", I3)') info
   STOP
ENDIF

! Generate the DVR grid, the weights and the KE matrix
ALLOCATE(gridpts(1:nbasis), gridweights(1:nbasis))
ALLOCATE(kinetic(1:nbasis,1:nbasis))
gridpts = czero
gridweights = czero
kinetic = czero

CALL dvr(numelements, celements, norder, nbasis, angle, gridpts, gridweights, kinetic)

! If call the insertion subroutine
! If kedvr = 1, we're using the separable representation of the kinentic energy
! If switchv > 1, we're using the separable representation of the potential energy Vnuc+2J-K
!       If switchv = 1 we're using the Pseudoinverse to orthongalize
!       If switchv = 2, we're using orthonormal orbitals to orthogonalize
! If partitionflag = 1, we're using a modified potential of vnuc by spliting it with P(r)
!  that is -Z/r = -Z/r * (1-P(r)) + -Z/r * P(r)
ALLOCATE(hamiltonian(1:nbasis,1:nbasis))
hamiltonian = czero
ALLOCATE(potential(1:nbasis,1:nbasis))
potential = czero

! If we're inserting the kinetic energy, we'll do it with the potential energy
IF(kedvr .EQ. 1) THEN
   kinetic = czero
ENDIF

!Generate the potential either on DVRs (switchv=0) or via insertion (switchv.NE.1)
IF(switchv .EQ. 0) THEN
   DO i=1,nbasis
      potential(i,i) = -1.0/gridpts(i) +0.5d0*(l_ang*(l_ang+1.0d0))/gridpts(i)/gridpts(i)
   ENDDO

! If insertion
ELSE
   CALL insertion(potential, gridpts, gridweights)

   ! Add the partitioning of the nuclear potential if partitionflag == 1
   IF(partitionflag .EQ. 1) THEN
      WRITE(*,*) "Using the partitioning function!"
      ALLOCATE(temppot(1:nbasis,1:nbasis))
      temppot = czero

      ! For the first element, use the exact nuclear potential for the first element
      DO i=1, norder-1
         temppot(i,i) = -10.0/gridpts(i)
      ENDDO

      ! For the second element, use the representation of -Z/r * (1-P(r))
      DO i=norder, 2*norder-2
         temppot(i,i) = -10.0*(1.0d0-partitioning(REAL(gridpts(i)), REAL(gridpts(norder-1)), REAL(gridpts(2*norder-2))))/gridpts(i)
      ENDDO      

      potential = potential + temppot
      DEALLOCATE(temppot)
   ELSE
      !Add the L^2 Eigenvalue
      IF (l_ang > 0) THEN
         DO i=1, nbasis
            potential(i,i) = potential(i,i) + 0.5d0*(l_ang*(l_ang+1.0d0))/gridpts(i)/gridpts(i)
         ENDDO
      ENDIF
   ENDIF
ENDIF

hamiltonian = kinetic + potential

!Calculates the energy spectrum
WRITE(*,*) " Calcuation ENERGY!"

ALLOCATE(wavefunction(1:nbasis))

CALL engspec(hamiltonian, wavefunction, wfn)

!!$! Test properties of the eigenvalues and eigenvectors
!!$CALL eigentest(kinetic, potential,gridpts,gridweights, wavefunction)


DEALLOCATE(gridpts, kinetic)
DEALLOCATE(gridweights)
DEALLOCATE(potential)
DEALLOCATE(wavefunction)
DEALLOCATE(hamiltonian)
ENDPROGRAM
