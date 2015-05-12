!This subroutine constructs the matrix of the Laguerre potential energy
!from the old code's lower diagonal matrix
SUBROUTINE svd_lagpot(overlaps, matrixin)
USE inputvariables
IMPLICIT NONE

INTEGER:: i,j
REAL*8 :: matrixin(1:numgauss, 1:numgauss)
REAL*8 :: kinetic(1:numgauss,1:numgauss)
REAL*8 :: direct(1:numgauss,1:numgauss)
REAL*8 :: exchange(1:numgauss,1:numgauss)
REAL*8 :: vnuc(1:numgauss,1:numgauss)
REAL*8 :: overlaps(1:numgauss, 1:numgauss)

!!$! DEBUG
!!$INTEGER :: lwork, info, ierror
!!$REAL(KIND=8) :: w(1:numgauss), work(1:10*numgauss)
!!$lwork=10*numgauss


! Read in overlaps
CALL readmesa('overlaps.dat', numgauss, overlaps)

! Read in kinetic energy
CALL readmesa('kineticx.dat', numgauss, kinetic)

! Read in vnuc
CALL readmesa('vnucxxxx.dat', numgauss, vnuc)

! Read in direct
CALL readmesa('directxx.dat', numgauss, direct)

! Read in exchange
CALL readmesa('exchange.dat', numgauss, exchange)

matrixin = vnuc + direct -0.5d0*exchange
!matrixin = vnuc 
!matrixin = vnuc + direct

!!$OPEN(UNIT=12, FILE="overlaps.out", STATUS='UNKNOWN', ACTION='WRITE')
!!$DO i=1, numgauss
!!$   WRITE(12,'(100ES19.9)') (overlaps(i,j),j=1, numgauss)
!!$ENDDO
!!$CLOSE(12)

!!$OPEN(UNIT=9, FILE='diagdebug.dat', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)
!!$CALL DSYGV(1, 'N', 'L', numgauss, matrixin, numgauss, overlaps, numgauss, w, work, lwork, info)
!!$DO i=1, numgauss
!!$   WRITE(9,'(ES17.9)') w(i)
!!$ENDDO
!!$CLOSE(9)
!!$
!!$STOP

ENDSUBROUTINE
