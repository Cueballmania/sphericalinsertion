!This subroutine evalues the integral, in the DVR approximation, of a DVR function and 
! an unnormalized gaussian
SUBROUTINE gaussmat(gridpts, gridweights, gaussdvr, mfac, part)
USE inputvariables
IMPLICIT NONE

COMPLEX(KIND=DBL), INTENT(IN) :: gridpts(1:nbasis)
COMPLEX(KIND=DBL), INTENT(IN) :: gridweights(1:nbasis)
INTEGER, INTENT(IN) :: mfac                                    ! m-factor for conjugation
INTEGER, INTENT(IN) :: part ! Partition space?

COMPLEX(KIND=DBL), INTENT(OUT) :: gaussdvr(1:nbasis,1:numprimg)  ! DVR matrix evaluation

INTEGER :: i, j, k
COMPLEX(KIND=DBL) :: dvrpt, dvrweight
COMPLEX(KIND=DBL) :: gaussarray(1:numprimg)
REAL(KIND=DBL) :: escale
REAL(KIND=DBL) :: offset
REAL(KIND=DBL) :: partition, expon

If (part .EQ. 0) THEN
   !Main loop
   DO i=1, nbasis
      dvrpt = gridpts(i)
      dvrweight = gridweights(i)

      CALL stos(dvrpt, numprimg, gaussarray, mfac, l_ang, m_ang)
      
      DO j=1, numprimg
         gaussdvr(i,j) = SQRT(dvrweight)*gaussarray(j)
      ENDDO
   ENDDO


ELSE
   DO i=1, norder-1
      DO j=1, numprimg
         gaussdvr(i,j) = czero
      ENDDO
   ENDDO

   escale = gridpts(2*norder-2) - gridpts(norder-1)
   offset = gridpts(norder-1) + escale
   
   DO i=norder, 2*norder-2
      dvrpt = gridpts(i)
      dvrweight = gridweights(i)
      expon = ((dvrpt-offset)/escale)
      partition = DEXP(expon*expon/(expon*expon-1))
      
      CALL stos(dvrpt, numprimg, gaussarray, mfac, l_ang, m_ang)
      
      DO j=1, numprimg
         gaussdvr(i,j) = SQRT(dvrweight)*gaussarray(j)*SQRT(partition)
      ENDDO
   ENDDO
   
   DO i=2*norder-1, nbasis
      dvrpt = gridpts(i)
      dvrweight = gridweights(i)

      CALL stos(dvrpt, numprimg, gaussarray, mfac, l_ang, m_ang)
      
      DO j=1, numprimg
         gaussdvr(i,j) = SQRT(dvrweight)*gaussarray(j)
      ENDDO
   ENDDO
ENDIF

!
!!$OPEN(unit=1000, file='gauss_dvr.out', status='unknown')
!!$DO i=1, nbasis
!!$   WRITE(1000,'(1x,100ES19.9E3)') DREAL(gridpts(i)), (gaussexp(j), lagmat(i,j),j=1,numgauss)
!!$ENDDO
!!$!
!!$CLOSE(1000)
!
END SUBROUTINE
