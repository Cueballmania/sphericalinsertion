!This subroutine makes an array of the potential
!The potentail dies off to zero on the complex scale
SUBROUTINE cdiagpot(cpotential, gridpts)
USE inputvariables
IMPLICIT NONE
!
INTEGER :: i, j
REAL(KIND=8) :: cpotential(1:nbasis)
COMPLEX(KIND=8) :: gridpts(1:nbasis)
!
OPEN(unit=1005, file='cenfpotent.out', status='unknown')
!
DO i=1, (norder-1)*relements
      cpotential(i)=DREAL(angular(angular+1))*5D-1/(DREAL(gridpts(i))*DREAL(gridpts(i))*mass)
      WRITE(1005,*) DREAL(gridpts(i)), cpotential(i)
ENDDO
!
DO i=(norder-1)*relements+1, nbasis
      cpotential(i)=0.0D0
      WRITE(1005,*) 0d0, cpotential(i)
ENDDO
!
CLOSE(1005)
!
ENDSUBROUTINE
