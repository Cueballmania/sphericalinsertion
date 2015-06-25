!This subroutine reads in the MESA matrices.
SUBROUTINE gausspot(numgauss,overlaps, matrixin, kedvr, partitionflag)
IMPLICIT NONE

! Double Precision
INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=13,r=200)

! Input variable
INTEGER, INTENT(IN) :: numgauss
INTEGER, INTENT(IN) :: kedvr
INTEGER, INTENT(IN) :: partitionflag

! Output matrices
REAL(KIND=DBL), INTENT(OUT) :: matrixin(1:numgauss, 1:numgauss)
REAL(KIND=DBL), INTENT(OUT) :: overlaps(1:numgauss, 1:numgauss)

! Local variables
INTEGER:: i,j
REAL(KIND=DBL) :: kinetic(1:numgauss,1:numgauss)
REAL(KIND=DBL) :: direct(1:numgauss,1:numgauss)
REAL(KIND=DBL) :: exchange(1:numgauss,1:numgauss)
REAL(KIND=DBL) :: vnuc(1:numgauss,1:numgauss)

!!!$DEBUG
!!$INTEGER :: lwork, info, ierror
!!$REAL(KIND=8) :: w(1:numgauss), work(1:10*numgauss)
!!$lwork=10*numgauss

! Read in overlaps
CALL readmesa('overlaps.dat', numgauss, overlaps)

! Read in kinetic energy
CALL readmesa('kineticx.dat', numgauss, kinetic)

! Read in vnuc
IF(partitionflag .EQ. 1) THEN
   CALL readmesa('vnucthet.dat', numgauss, vnuc)
ELSEIF(partitionflag .EQ. 2) THEN
   CALL readmesa('vnucbeck.dat', numgauss, vnuc)
ELSEIF(partitionflag .EQ. 3) THEN
   CALL readmesa('vnucbump.dat', numgauss, vnuc) 
ELSE
   CALL readmesa('vnucxxxx.dat', numgauss, vnuc)
ENDIF

! Read in direct
CALL readmesa('directxx.dat', numgauss, direct)

! Read in exchange
CALL readmesa('exchange.dat', numgauss, exchange)

! Include vnuc?
IF (kedvr .EQ. 1) THEN
   matrixin =  kinetic + vnuc + direct - 0.5*exchange
ELSE
   matrixin = vnuc + direct - 0.5*exchange
ENDIF


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
