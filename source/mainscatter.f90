MODULE inputvariables
IMPLICIT NONE
SAVE

! Declare double precision
INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=13,r=200)

! Maximum number of quadrature points per element
INTEGER, PARAMETER :: MAXN = 500

! Static parameters
!-- Nuclear Charge could become an input parameter -- !
REAL(KIND=DBL), PARAMETER :: PI = 3.141592653589
INTEGER, PARAMETER :: nucharge = 4
REAL(KIND=DBL), PARAMETER ::TOLERANCE = 1D-10
REAL(KIND=DBL), PARAMETER :: duey=2.0d0
COMPLEX(KIND=DBL), PARAMETER :: czero=(0.0D0,0.0D0)
COMPLEX(KIND=DBL), PARAMETER :: cone=(1.0D0,0.0D0)
COMPLEX(KIND=DBL), PARAMETER :: cimag=(0.0D0,1.0D0)
REAL(KIND=DBL), PARAMETER :: mass=1.0

! Quadrature gobals
INTEGER :: numelements, celements, relements  ! Number of elements, number of complex elements, number of real elements
INTEGER :: norder, angular

! DVR gobals
INTEGER :: nbasis
REAL(KIND=DBL) :: gridend
COMPLEX(KIND=DBL) :: contourphase       ! Complex number of the real scaling

! Insertion gobals
INTEGER :: numgauss, numprimg

! Calculation inputs
INTEGER :: kedvr                ! Determines how the Kinetic Energy is used
                                ! 0: DVR generated KE (with centifugal term if needed)
                                ! 1: Input from kineticx.dat generated from MESA in terms of the contracted Gaussians

INTEGER :: switchv              ! Determines how the Potential is used
                                ! 0: DVR generated potential that is the bare nuclear attraction potential (solves the hydrogen-like problem)
                                ! 1: Uses insertion -- possibly with partitioning the nuclear potential(depends on partitionflag) -- with
                                !         the pseudoinverse of the SVD.
                                ! 2: Same as 1 except uses the orthogonalization of the SVD orbitals -- produces the same results as 1.
                                ! 3: Reads in only the Exchange matrix from the mesa calculation.  The coefficents of the contracted Gaussians
                                !          for occupied orbitals must be hardcoded in direct_coul.f90.  The nuclear potential is generated on the DVR

INTEGER :: l_ang, m_ang         ! Angular values

INTEGER :: partitionflag        ! Partitioning flag to partition the nuclear attraction potential
                                ! 0: no paritioning
                                ! 1: Heaviside-theta paritioning. Zero in the first two elements.  Reads from vnucthet.dat
                                ! 2: 3rd order Becke paritioning in the second element. Reads from vnucbeck.dat
                                ! 3: Bump function partitioning in the second element.  Reads from vnucbump.dat

REAL(KIND=DBL) :: SVD_tol       ! Min value for SVD
REAL(KIND=DBL) :: contourpt     ! Point where scaling starts
REAL(KIND=DBL) :: angle         ! Angle in degrees
END MODULE

! This program was part of the PhD work of Brant Abeln.  The bulk of this code was coded by him, but he was not without help.
! Many people collaborated and gave ideas towards the completion of his project.
! This code represents the spherical polar coordinate version of the separable insertion method.
! It aims to find static-exchange resonance eigenfunction and eigenvalues for atoms while also recalculating the occupied orbital energies and functions.
!


PROGRAM mainscatter
USE inputvariables
IMPLICIT NONE

! Local variables
INTEGER :: i, j, info

! Hamiltonian matrices
COMPLEX(KIND=DBL), ALLOCATABLE :: kinetic(:,:), potential(:,:)
COMPLEX(KIND=DBL), ALLOCATABLE ::  temppot(:,:), temppot2(:)
COMPLEX(KIND=DBL), ALLOCATABLE :: hamiltonian(:,:)

! Basis information
COMPLEX(KIND=DBL), ALLOCATABLE :: gridpts(:)
COMPLEX(KIND=DBL), ALLOCATABLE :: gridweights(:)

! Wavefunction test of properites of the ground state
INTEGER :: wfn
COMPLEX(KIND=DBL), ALLOCATABLE :: wavefunction(:)

! Partition function
REAL(KIND=DBL) :: beckepart, bumppart
REAL(KIND=DBL) :: parfact

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
nbasis = (norder-1)*numelements
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

! Add the centrifugal term
IF (l_ang > 0 .AND. (kedvr .NE. 1 )) THEN
   DO i=1, nbasis
      kinetic(i,i) = kinetic(i,i) + 0.5d0*(l_ang*(l_ang+1.0d0))/gridpts(i)/gridpts(i)
   ENDDO
ENDIF

!Generate the potential either on DVRs (switchv=0) or via insertion (switchv.NE.1)
IF(switchv .EQ. 0) THEN
   DO i=1,nbasis
      potential(i,i) = -1.0/gridpts(i) 
   ENDDO

! If insertion, check if we're making a partitioned potential
ELSE IF(switchv .EQ. 1 .OR. switchv .EQ. 2) THEN
   ! Add the partitioning of the nuclear potential if partitionflag == 1
   IF(partitionflag .GT. 0) THEN
      WRITE(*,*) "Using the partitioning function!"
      WRITE(*,*) "Creating the partitioned nuclear attraction matrix"

!      CALL vnucpart(nucharge, norder, numelements, numgauss, numprimg)

      WRITE(*,*) "Finished creating partitioned nuclear attraction matrix"

      CALL insertion(potential, gridpts, gridweights)
      
      ALLOCATE(temppot(1:nbasis,1:nbasis))
      temppot = czero

      ! For the first element, use the exact nuclear potential for the first element
      DO i=1, norder-1
         temppot(i,i) = -nucharge/gridpts(i)
      ENDDO

      ! For the second element, use the representation of -Z/r * (1-P(r))
      DO i=norder, 2*norder-2
         ! Heviside Theta function
         IF(partitionflag .EQ. 1) THEN
            parfact = 0.0d0
         ! Becke partitioning function
         ELSEIF(partitionflag .EQ. 2) THEN
            parfact = (1.0d0-beckepart(REAL(gridpts(i)), REAL(gridpts(norder-1)), REAL(gridpts(2*norder-1))))
         ! Bump partitioning function
         ELSEIF(partitionflag .EQ. 3) THEN
            parfact = (1.0d0-bumppart(REAL(gridpts(i)), REAL(gridpts(norder-1)), REAL(gridpts(2*norder-1))))
         ELSE
            WRITE(*,*) "parititionflag = ", partitionflag, " not defined"
         ENDIF
         WRITE(888,*) parfact, REAL(gridpts(i))
         temppot(i,i) = -nucharge*parfact/gridpts(i)
      ENDDO

      OPEN(UNIT=23, FILE='potentialcut.out', STATUS='UNKNOWN', ACTION='WRITE')
      DO i=1, nbasis
         WRITE(23, *) REAL(gridpts(i)), IMAG(gridpts(i)), REAL(temppot(i,i))
      ENDDO
      CLOSE(23)
      
      potential = potential + temppot
      DEALLOCATE(temppot)

   ! If using regular insertion for the nuclear attraction term just do it.
   ELSE
      CALL insertion(potential, gridpts, gridweights)
   ENDIF

! Directly calculate the 2J term, the Vnuc term and insertion with the K matrix
! The 2J term is hard-coded in direct_coul.f90   
ELSE IF(switchv .EQ. 3) THEN
   WRITE(*,*) "Separable-Exchange calculation"
   WRITE(*,*) "Putting the direct Coulomb potential directly on the grid"
   WRITE(*,*) "****Assuming Beryllium -- other systems not implemented yet!!!"

   ALLOCATE(temppot(1:nbasis,1:nbasis),temppot2(1:nbasis))
   temppot = czero
   temppot2 = czero

   ! Use the exact nuclear potential for the DVR
   DO i=1, nbasis
      temppot(i,i) = -nucharge/gridpts(i)
   ENDDO

   OPEN(UNIT=666, file='nucatt.out', status='unknown')
   DO i=1, nbasis
      WRITE(666,'(1x,4ES17.9)') DREAL(gridpts(i)), DIMAG(gridpts(i)), DREAL(temppot(i,i)), DIMAG(temppot(i,i))
   ENDDO
   CLOSE(666)
   
   ! Generate the direct Coulomb potential
   CALL direct_coul(gridpts,temppot2)

   ! Insertion where hopefully the exchange matrix is the only nonzero matrix.
   CALL insertion(potential, gridpts, gridweights)

   DO i=1, nbasis
      potential(i,i) = potential(i,i) + temppot(i,i) + 2.0*temppot2(i)
   ENDDO

ENDIF

! Form the Hamiltonian matrix
hamiltonian = kinetic + potential

!Calculates the energy spectrum
WRITE(*,*) " Calcuation ENERGY!"

ALLOCATE(wavefunction(1:nbasis))

! Calculate the spectrum and kick back the wfn-th eigenvalue
CALL engspec(hamiltonian, wavefunction, wfn)

IF(wfn .NE. 0)THEN
   OPEN(UNIT=222, file='wavefunction.out', status='unknown')
   DO i=1, nbasis
      WRITE(222,'(1x,4ES17.9)') DREAL(gridpts(i)), DIMAG(gridpts(i)), &
      & DREAL(wavefunction(i)/SQRT(gridweights(i))), DIMAG(wavefunction(i)/SQRT(gridweights(i)))
   ENDDO
   CLOSE(222)
ENDIF

!!$! Test properties of the eigenvalues and eigenvectors
!!$CALL eigentest(kinetic, potential,gridpts,gridweights, wavefunction)


DEALLOCATE(gridpts, kinetic)
DEALLOCATE(gridweights)
DEALLOCATE(potential)
DEALLOCATE(wavefunction)
DEALLOCATE(hamiltonian)
ENDPROGRAM
