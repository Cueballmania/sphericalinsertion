SUBROUTINE stos(gridpt, nsize, garray, mfac, lang, mang)
! This subroutine evalues the value of all of the primitive
! gaussian functions at a particular point and multiples
! the exact angular integral
! Note:
! x = r sin(theta) cos(phi)
! y = r sin(theta) sin(phi)
! z = r cos(theta)
IMPLICIT NONE

! Double Precision
INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=13,r=200)

! Useful angular constants
REAL(KIND=DBL), PARAMETER :: pi = 3.14159265358979
REAL(KIND=DBL),PARAMETER :: sqtpi = 1.77245385090552
COMPLEX(KIND=DBL), PARAMETER :: ci = (0.0d0,1.0d0)

! Input variables
INTEGER, INTENT(IN) :: nsize
INTEGER, INTENT(IN) :: mfac                 ! if mfac=1, Y; if mfac=-1, Y*
INTEGER, INTENT(IN) :: lang, mang
COMPLEX(KIND=DBL), INTENT(IN) :: gridpt

! Output array
COMPLEX(KIND=DBL), INTENT(OUT) :: garray(1:nsize)

! Local variables
INTEGER :: i,j, iread
INTEGER :: ierror
REAL(KIND=DBL) :: expo
COMPLEX(KIND=DBL):: evalue
CHARACTER(LEN=2) :: sym
REAL(KIND=DBL) :: zoff

! Open the file!
OPEN(UNIT=10, FILE="expos.dat", STATUS="OLD", ACTION="READ", IOSTAT=ierror)
fileopen: IF (ierror == 0) THEN

   readeval:DO i=1, nsize
      
      ! Read the ith entry in expos.dat
      READ(10,*,IOSTAT=iread) zoff, sym, expo

      ! Print message if there is a read error
      IF (iread /= 0) THEN
         WRITE(*,100) i, REAL(gridpt), AIMAG(gridpt)
100      FORMAT (1x,"Error reading entry: ", I3, " from expos.dat for r= ", ES15.6, " +i ", ES15.6)
         WRITE(*,101) zoff, sym, expo, iread
101      FORMAT (1x,"zoff= ", ES15.6, " sym= ", A, " expo= ", ES15.6, " iread= ", I10)
      ENDIF
      
      ! Evaluate the Gaussian expoential
      evalue = EXP(-expo*(gridpt*gridpt))

      ! Cases for the angular parts
      IF(sym .EQ. 'ss' .AND. lang .EQ. 0) THEN
         evalue = (2.0d0*expo/pi)**(0.75)*gridpt*evalue*2.0d0*sqtpi

      ELSEIF(sym .EQ. 'px' .AND. lang .EQ. 1 .AND. mang .EQ. 1) THEN
         evalue = -gridpt*gridpt*evalue*SQRT(2.0d0/3.0d0)*sqtpi*2.0d0*(2.0d0/pi)**(0.75)*expo**(1.25)
      ELSEIF(sym .EQ. 'px' .AND. lang .EQ. 1 .AND. mang .EQ. -1) THEN
         evalue = gridpt*gridpt*evalue*SQRT(2.0d0/3.0d0)*sqtpi*2.0d0*(2.0d0/pi)**(0.75)*expo**(1.25)
      ELSEIF(sym .EQ. 'py' .AND. lang .EQ. 1 .AND. ABS(mang) .EQ. 1 .AND. mfac .EQ. 1) THEN
         evalue = -ci*gridpt*gridpt*evalue*SQRT(2.0d0/3.0d0)*sqtpi*2.0d0*(2.0d0/pi)**(0.75)*expo**(1.25)
      ELSEIF(sym .EQ. 'py' .AND. lang .EQ. 1 .AND. ABS(mang) .EQ. 1 .AND. mfac .EQ. -1) THEN
         evalue = ci*gridpt*gridpt*evalue*SQRT(2.0d0/3.0d0)*sqtpi*2.0d0*(2.0d0/pi)**(0.75)*expo**(1.25)
      ELSEIF(sym .EQ. 'pz' .AND. lang .EQ. 1 .AND. mang .EQ. 0) THEN
         evalue = gridpt*gridpt*evalue*2.0d0*sqtpi/SQRT(3.0d0)*2.0d0*(2.0d0/pi)**(0.75)*expo**(1.25)

      ELSE
         evalue = (0.0d0,0.0d0)
      ENDIF

      garray(i) = evalue
      WRITE(1909,'(1x,I3,3ES16.7E2,1x,A3,5(ES15.6E3,1x))') mang,gridpt,zoff,sym,expo,evalue
   ENDDO readeval

CLOSE(10)
! If there was a read error, let me know!
ELSE fileopen
   WRITE(*,*) " Error opening expos.dat"
ENDIF fileopen
END SUBROUTINE
