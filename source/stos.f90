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
CHARACTER(LEN=3) :: sym
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
      evalue = (2.0d0*expo/pi)**(0.75)*EXP(-expo*(gridpt*gridpt))

      ! Cases for the angular parts
      IF(sym .EQ. 'ss' .AND. lang .EQ. 0) THEN
         evalue = evalue*2.0d0*sqtpi

      ELSEIF(sym .EQ. 'px' .AND. lang .EQ. 1 .AND. mang .EQ. 1) THEN
         evalue = -gridpt*evalue*SQRT(2.0d0/3.0d0)*sqtpi*2.0d0*SQRT(expo)
      ELSEIF(sym .EQ. 'px' .AND. lang .EQ. 1 .AND. mang .EQ. -1) THEN
         evalue = gridpt*evalue*SQRT(2.0d0/3.0d0)*sqtpi*2.0d0*SQRT(expo)
      ELSEIF(sym .EQ. 'py' .AND. lang .EQ. 1 .AND. ABS(mang) .EQ. 1 .AND. mfac .EQ. 1) THEN
         evalue = -ci*gridpt*evalue*SQRT(2.0d0/3.0d0)*sqtpi*2.0d0*SQRT(expo)
      ELSEIF(sym .EQ. 'py' .AND. lang .EQ. 1 .AND. ABS(mang) .EQ. 1 .AND. mfac .EQ. -1) THEN
         evalue = ci*gridpt*evalue*SQRT(2.0d0/3.0d0)*sqtpi*2.0d0*SQRT(expo)
      ELSEIF(sym .EQ. 'pz' .AND. lang .EQ. 1 .AND. mang .EQ. 0) THEN
         evalue = gridpt*evalue*4.0d0*sqtpi*SQRT(expo/3.0d0)


      ! Normalization of dxx, dyy, dzz, functions are off by a factor of sqrt(3)
      ELSEIF(sym .EQ. 'xx' .AND. lang .EQ. 0 .AND. mang .EQ. 0) THEN
         evalue = gridpt*gridpt*8.0d0*expo*sqtpi/3.0d0*evalue
      ELSEIF(sym .EQ. 'xx' .AND. lang .EQ. 2 .AND. mang .EQ. 0) THEN
         evalue = -gridpt*gridpt*8.0d0*expo*sqtpi/3.0d0/SQRT(5.0d0)*evalue
      ELSEIF(sym .EQ. 'xx' .AND. lang .EQ. 2 .AND. ABS(mang) .EQ. 2) THEN
         evalue = gridpt*gridpt*4.0d0*expo*sqtpi*SQRT(2.0d0/15.0d0)*evalue
      ELSEIF(sym .EQ. 'yy' .AND. lang .EQ. 0 .AND. mang .EQ. 0) THEN
         evalue = gridpt*gridpt*8.0d0*expo*sqtpi/3.0d0*evalue
      ELSEIF(sym .EQ. 'yy' .AND. lang .EQ. 2 .AND. mang .EQ. 0) THEN
         evalue = -gridpt*gridpt*8.0d0*expo*sqtpi/3.0d0/SQRT(5.0d0)*evalue
      ELSEIF(sym .EQ. 'yy' .AND. lang .EQ. 2 .AND. ABS(mang) .EQ. 2) THEN
         evalue = -gridpt*gridpt*4.0d0*expo*sqtpi*SQRT(2.0d0/15.0d0)*evalue
      ELSEIF(sym .EQ. 'zz' .AND. lang .EQ. 0 .AND. mang .EQ. 0) THEN
         evalue = gridpt*gridpt*8.0d0*expo*sqtpi/3.0d0*evalue
      ELSEIF(sym .EQ. 'zz' .AND. lang .EQ. 2 .AND. mang .EQ. 0) THEN
         evalue = gridpt*gridpt*16.0d0*expo*sqtpi/3.0d0/SQRT(5.0d0)*evalue

      ELSEIF(sym .EQ. 'xy' .AND. lang .EQ. 2 .AND. mang .EQ. 2 .AND. mfac .EQ. 1) THEN
         evalue = ci*gridpt*gridpt*4.0d0*expo*SQRT(2.0d0/15.0d0)*sqtpi*evalue
      ELSEIF(sym .EQ. 'xy' .AND. lang .EQ. 2 .AND. mang .EQ. 2 .AND. mfac .EQ. -1) THEN
         evalue = -ci*gridpt*gridpt*4.0d0*expo*SQRT(2.0d0/15.0d0)*sqtpi*evalue
      ELSEIF(sym .EQ. 'xy' .AND. lang .EQ. 2 .AND. mang .EQ. -2 .AND. mfac .EQ. 1) THEN
         evalue = -ci*gridpt*gridpt*4.0d0*expo*SQRT(2.0d0/15.0d0)*sqtpi*evalue
      ELSEIF(sym .EQ. 'xy' .AND. lang .EQ. 2 .AND. mang .EQ. -2 .AND. mfac .EQ. -1) THEN
         evalue = ci*gridpt*gridpt*4.0d0*expo*SQRT(2.0d0/15.0d0)*sqtpi*evalue         
      ELSEIF(sym .EQ. 'xz' .AND. lang .EQ. 2 .AND. mang .EQ. 1) THEN
         evalue = -gridpt*gridpt*4.0d0*expo*SQRT(2.0d0/15.0d0)*sqtpi*evalue
      ELSEIF(sym .EQ. 'xz' .AND. lang .EQ. 2 .AND. mang .EQ. -1) THEN
         evalue = gridpt*gridpt*4.0d0*expo*SQRT(2.0d0/15.0d0)*sqtpi*evalue
      ELSEIF(sym .EQ. 'yz' .AND. lang .EQ. 2 .AND. ABS(mang) .EQ. 1 .AND. mfac .EQ. 1) THEN
         evalue = -ci*gridpt*gridpt*4.0d0*expo*SQRT(2.0d0/15.0d0)*sqtpi*evalue
      ELSEIF(sym .EQ. 'yz' .AND. lang .EQ. 2 .AND. ABS(mang) .EQ. 1 .AND. mfac .EQ. -1) THEN
         evalue = ci*gridpt*gridpt*4.0d0*expo*SQRT(2.0d0/15.0d0)*sqtpi*evalue


      ! f-functions normalized for xyz; off by sqrt(3) for nnm and off by sqrt(15) for nnn functions
      ELSEIF(sym .EQ. 'xxx' .AND. lang .EQ. 3 .AND. mang .EQ. 3) THEN
         evalue = -gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi/SQRT(35.0d0)*evalue
      ELSEIF(sym .EQ. 'xxx' .AND. lang .EQ. 3 .AND. mang .EQ. -3) THEN
         evalue = gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi/SQRT(35.0d0)*evalue
      ELSEIF(sym .EQ. 'xxx' .AND. lang .EQ. 3 .AND. mang .EQ. 1) THEN
         evalue = gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi*SQRT(3.0d0/7.0d0)*0.2d0*evalue
      ELSEIF(sym .EQ. 'xxx' .AND. lang .EQ. 3 .AND. mang .EQ. -1) THEN
         evalue = -gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi*SQRT(3.0d0/7.0d0)*0.2d0*evalue
      ELSEIF(sym .EQ. 'xxx' .AND. lang .EQ. 1 .AND. mang .EQ. 1) THEN
         evalue = -gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi*SQRT(6.0d0)*0.2d0*evalue
      ELSEIF(sym .EQ. 'xxx' .AND. lang .EQ. 1 .AND. mang .EQ. -1) THEN
         evalue = gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi*SQRT(6.0d0)*0.2d0*evalue

      ELSEIF(sym .EQ. 'yyy' .AND. lang .EQ. 3 .AND. mang .EQ. 3 .AND. mfac .EQ. 1) THEN
         evalue = -ci*gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi/SQRT(35.0d0)*evalue
      ELSEIF(sym .EQ. 'yyy' .AND. lang .EQ. 3 .AND. mang .EQ. 3 .AND. mfac .EQ. -1) THEN
         evalue = ci*gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi/SQRT(35.0d0)*evalue
      ELSEIF(sym .EQ. 'yyy' .AND. lang .EQ. 3 .AND. mang .EQ. -3 .AND. mfac .EQ. 1) THEN
         evalue = ci*gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi/SQRT(35.0d0)*evalue
      ELSEIF(sym .EQ. 'yyy' .AND. lang .EQ. 3 .AND. mang .EQ. -3 .AND. mfac .EQ. -1) THEN
         evalue = -ci*gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi/SQRT(35.0d0)*evalue
      ELSEIF(sym .EQ. 'yyy' .AND. lang .EQ. 3 .AND. ABS(mang) .EQ. 1 .AND. mfac .EQ. 1) THEN
         evalue = ci*gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi*SQRT(3.0d0/7.0d0)*0.2d0*evalue
      ELSEIF(sym .EQ. 'yyy' .AND. lang .EQ. 3 .AND. ABS(mang) .EQ. 1 .AND. mfac .EQ. -1) THEN
         evalue = -ci*gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi*SQRT(3.0d0/7.0d0)*0.2d0*evalue
      ELSEIF(sym .EQ. 'yyy' .AND. lang .EQ. 1 .AND. ABS(mang) .EQ. 1 .AND. mfac .EQ. 1) THEN
         evalue = -ci*gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi*SQRT(6.0d0)*0.2d0*evalue
      ELSEIF(sym .EQ. 'yyy' .AND. lang .EQ. 1 .AND. ABS(mang) .EQ. 1 .AND. mfac .EQ. -1) THEN
         evalue = ci*gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi*SQRT(6.0d0)*0.2d0*evalue

      ELSEIF(sym .EQ. 'zzz' .AND. lang .EQ. 3 .AND. mang .EQ. 0) THEN
         evalue = gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi/SQRT(7.0d0)*0.8d0*evalue
      ELSEIF(sym .EQ. 'zzz' .AND. lang .EQ. 1 .AND. mang .EQ. 0) THEN
         evalue = gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi*SQRT(3.0d0)*0.4d0*evalue

      ELSEIF(sym .EQ. 'xxy' .AND. lang .EQ. 3 .AND. ABS(mang) .EQ. 3 .AND. mfac .EQ. 1) THEN
         evalue = -ci*gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi/SQRT(35.0d0)*evalue
      ELSEIF(sym .EQ. 'xxy' .AND. lang .EQ. 3 .AND. ABS(mang) .EQ. 3 .AND. mfac .EQ. -1) THEN
         evalue = ci*gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi/SQRT(35.0d0)*evalue
      ELSEIF(sym .EQ. 'xxy' .AND. lang .EQ. 3 .AND. ABS(mang) .EQ. 1 .AND. mfac .EQ. 1) THEN
         evalue = ci*gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi/SQRT(21.0d0)*0.20d0*evalue
      ELSEIF(sym .EQ. 'xxy' .AND. lang .EQ. 3 .AND. ABS(mang) .EQ. 1 .AND. mfac .EQ. -1) THEN
         evalue = -ci*gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi/SQRT(21.0d0)*0.20d0*evalue
      ELSEIF(sym .EQ. 'xxy' .AND. lang .EQ. 1 .AND. ABS(mang) .EQ. 1 .AND. mfac .EQ. 1) THEN
         evalue = -ci*gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi*SQRT(2.0d0/3.0d0)*0.20d0*evalue
      ELSEIF(sym .EQ. 'xxy' .AND. lang .EQ. 1 .AND. ABS(mang) .EQ. 1 .AND. mfac .EQ. -1) THEN
         evalue = ci*gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi*SQRT(2.0d0/3.0d0)*0.20d0*evalue

      ELSEIF(sym .EQ. 'xxz' .AND. lang .EQ. 3 .AND. ABS(mang) .EQ. 2) THEN
         evalue = gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi*SQRT(2.0d0/105.0d0)*evalue
      ELSEIF(sym .EQ. 'xxz' .AND. lang .EQ. 3 .AND. mang .EQ. 0) THEN
         evalue = -gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi/SQRT(7.0d0)*0.4d0*evalue
      ELSEIF(sym .EQ. 'xxz' .AND. lang .EQ. 1 .AND. mang .EQ. 0) THEN
         evalue = gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi/SQRT(3.0d0)*0.4d0*evalue

      ELSEIF(sym .EQ. 'xyy' .AND. lang .EQ. 3 .AND. mang .EQ. 3) THEN
         evalue = gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi/SQRT(35.0d0)*evalue         
      ELSEIF(sym .EQ. 'xyy' .AND. lang .EQ. 3 .AND. mang .EQ. -3) THEN
         evalue = -gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi/SQRT(35.0d0)*evalue
      ELSEIF(sym .EQ. 'xyy' .AND. lang .EQ. 3 .AND. mang .EQ. 1) THEN
         evalue = gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi/SQRT(21.0d0)*0.2d0*evalue         
      ELSEIF(sym .EQ. 'xyy' .AND. lang .EQ. 3 .AND. mang .EQ. -1) THEN
         evalue = -gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi/SQRT(21.0d0)*0.2d0*evalue      
      ELSEIF(sym .EQ. 'xyy' .AND. lang .EQ. 1 .AND. mang .EQ. -1) THEN
         evalue = gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi*SQRT(2.0d0/3.0d0)*0.2d0*evalue         
      ELSEIF(sym .EQ. 'xyy' .AND. lang .EQ. 1 .AND. mang .EQ. 1) THEN
         evalue = -gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi*SQRT(2.0d0/3.0d0)*0.2d0*evalue       

      ELSEIF(sym .EQ. 'yyz' .AND. lang .EQ. 3 .AND. mang .EQ. 2) THEN
         evalue = gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi*SQRT(2.0d0/105.0d0)*evalue         
      ELSEIF(sym .EQ. 'yyz' .AND. lang .EQ. 3 .AND. mang .EQ. -2) THEN
         evalue = -gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi*SQRT(2.0d0/105.0d0)*evalue
      ELSEIF(sym .EQ. 'yyz' .AND. lang .EQ. 3 .AND. mang .EQ. 0) THEN
         evalue = -gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi/SQRT(7.0d0)*0.4d0*evalue
      ELSEIF(sym .EQ. 'yyz' .AND. lang .EQ. 1 .AND. mang .EQ. 0) THEN
         evalue = gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi/SQRT(3.0d0)*0.4d0*evalue

      ELSEIF(sym .EQ. 'xzz' .AND. lang .EQ. 3 .AND. mang .EQ. 1) THEN
         evalue = -gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi/SQRT(21.0d0)*0.8d0*evalue         
      ELSEIF(sym .EQ. 'xzz' .AND. lang .EQ. 3 .AND. mang .EQ. -1) THEN
         evalue = gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi/SQRT(21.0d0)*0.8d0*evalue
      ELSEIF(sym .EQ. 'xzz' .AND. lang .EQ. 1 .AND. mang .EQ. 1) THEN
         evalue = -gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi*SQRT(2.0d0/3.0d0)*0.2d0*evalue         
      ELSEIF(sym .EQ. 'xzz' .AND. lang .EQ. 1 .AND. mang .EQ. -1) THEN
         evalue = gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi*SQRT(2.0d0/3.0d0)*0.2d0*evalue   

      ELSEIF(sym .EQ. 'yzz' .AND. lang .EQ. 3 .AND. ABS(mang) .EQ. 1 .AND. mfac .EQ. 1) THEN
         evalue = -ci*gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi/SQRT(21.0d0)*0.8d0*evalue 
      ELSEIF(sym .EQ. 'yzz' .AND. lang .EQ. 3 .AND. ABS(mang) .EQ. 1 .AND. mfac .EQ. -1) THEN
         evalue = ci*gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi/SQRT(21.0d0)*0.8d0*evalue 
      ELSEIF(sym .EQ. 'yzz' .AND. lang .EQ. 1 .AND. ABS(mang) .EQ. 1 .AND. mfac .EQ. 1) THEN
         evalue = -ci*gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi*SQRT(2.0d0/3.0d0)*0.2d0*evalue 
      ELSEIF(sym .EQ. 'yzz' .AND. lang .EQ. 1 .AND. ABS(mang) .EQ. 1 .AND. mfac .EQ. -1) THEN
         evalue = ci*gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi*SQRT(2.0d0/3.0d0)*0.2d0*evalue

      ELSEIF(sym .EQ. 'xyz' .AND. lang .EQ. 3 .AND. mang .EQ. 2 .AND. mfac .EQ. 1) THEN
         evalue = ci*gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi*SQRT(2.0d0/105.0d0)*evalue
      ELSEIF(sym .EQ. 'xyz' .AND. lang .EQ. 3 .AND. mang .EQ. 2 .AND. mfac .EQ. -1) THEN
         evalue = -ci*gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi*SQRT(2.0d0/105.0d0)*evalue
      ELSEIF(sym .EQ. 'xyz' .AND. lang .EQ. 3 .AND. mang .EQ. -2 .AND. mfac .EQ. 1) THEN
         evalue = -ci*gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi*SQRT(2.0d0/105.0d0)*evalue
      ELSEIF(sym .EQ. 'xyz' .AND. lang .EQ. 3 .AND. mang .EQ. -2 .AND. mfac .EQ. -1) THEN
         evalue = ci*gridpt*gridpt*gridpt*8.0d0*expo**(1.5)*sqtpi*SQRT(2.0d0/105.0d0)*evalue
         

         
      ELSE
         evalue = (0.0d0,0.0d0)
      ENDIF

      garray(i) = evalue
      !WRITE(1909,'(1x,I3,3ES16.7E2,1x,A3,5(ES15.6E3,1x))') mang,gridpt,zoff,sym,expo,evalue
   ENDDO readeval

CLOSE(10)
! If there was a read error, let me know!
ELSE fileopen
   WRITE(*,*) " Error opening expos.dat"
ENDIF fileopen
END SUBROUTINE
