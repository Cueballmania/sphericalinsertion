!This subroutine will output the energy spectrum
SUBROUTINE engspec(hamiltonian, wavefunction, wfn)
USE inputvariables
IMPLICIT NONE

! Input variables
COMPLEX(KIND=DBL), INTENT(IN) :: hamiltonian(1:nbasis,1:nbasis)
INTEGER, INTENT(IN) :: wfn                                      ! eigenvector to test

! Output wavefuncton
COMPLEX(KIND=DBL), INTENT(OUT) :: wavefunction(1:nbasis)

! Local variables
INTEGER :: i, j

! LAPACK variables
INTEGER :: info
CHARACTER(LEN=1) :: jobvl = 'N', jobvr = 'V'
COMPLEX(KIND=DBL) :: energies(1:nbasis)
COMPLEX(KIND=DBL) :: empt(1:nbasis,1:nbasis)
COMPLEX(KIND=DBL) :: eigen(1:nbasis,1:nbasis)
COMPLEX(KIND=DBL) :: workspace(1:20*nbasis)
REAL(KIND=DBL) :: cwork(1:2*nbasis)


! Generate the hamiltonian matrix

!Diagonalize the hamiltonian matrix
CALL ZGEEV(jobvl, jobvr, nbasis, hamiltonian, nbasis, energies, empt, nbasis, eigen, nbasis, workspace, 20*nbasis, cwork, info)

WRITE(*,*) "Diagonal"

DO i=1, nbasis
   wavefunction(i) = eigen(i,wfn)
ENDDO

!Write out the energy spectrum
OPEN(unit=1000, file='spectrum.out', status='unknown')
DO i=1, nbasis
   WRITE(1000,'(1x,80ES17.9)') DREAL(energies(i)), DIMAG(energies(i))
ENDDO
CLOSE(1000)

!WRITE out the bound state spectrum
OPEN(unit=1001, file='bound.out', status='unknown')
DO i=1, nbasis
   IF(DREAL(energies(i)) .LE. -0.00d0) THEN
      WRITE(1001,'(1x,80ES17.9)') DREAL(energies(i)), DIMAG(energies(i))
      WRITE(*,'(1x,2ES17.9,1x,I4)') DREAL(energies(i)), DIMAG(energies(i)), i
   ENDIF
ENDDO
CLOSE(1001)

!!$!WRITE out the specturm in the range [3,4]
!!$OPEN(unit=1002, file='resrange.out', status='unknown')
!!$WRITE(*,*) "Res range"
!!$DO i=1, nbasis
!!$   IF(DREAL(energies(i)) .GE. 3.0d0 .AND. DREAL(energies(i)) .LE. 4.0d0) THEN
!!$      WRITE(1002,'(1x,10ES17.9)') DREAL(energies(i)), DIMAG(energies(i))
!!$      WRITE(*,'(1x,10ES17.9)') DREAL(energies(i)), DIMAG(energies(i))
!!$   ENDIF
!!$ENDDO
!!$CLOSE(1002)

ENDSUBROUTINE
