SUBROUTINE vnucpart(nucharge, norder, nelements, ngauss, nprim)
! This subroutine utilizes the underlying quadrature
! to re-calculate the nuclear attraction
! matrix elements using the partitioning function.
! This subroutine returns nothing, but generates the file
! vnucpart.dat
!
! This routine also reads FEM.in and xform.dat
USE gaussquad
IMPLICIT NONE

! Declare double precision
INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=13,r=200)
REAL(KIND=DBL), PARAMETER :: PI = 3.141592653589
REAL(KIND=DBL), PARAMETER :: tol = 1.0D-10

! Input variables
INTEGER, INTENT(IN) :: nucharge    ! >0
INTEGER, INTENT(IN) :: norder
INTEGER, INTENT(IN) :: nelements
INTEGER, INTENT(IN) :: ngauss
INTEGER, INTENT(IN) :: nprim

! Quadrature arrays
REAL(KIND=DBL), DIMENSION(1:norder) :: gausspts, gausswts, work
REAL(KIND=DBL):: elementsize(1:nelements)
REAL(KIND=DBL), ALLOCATABLE :: pts(:), dvrpts(:)
REAL(KIND=DBL), ALLOCATABLE :: wts(:), dvrwts(:)

! Gauss basis info
REAL(KIND=DBL) :: expos(1:nprim)
CHARACTER(LEN=2) :: sym(1:nprim)
REAL(KIND=DBL) :: xform(1:ngauss,1:nprim)

! Matrices
REAL(KIND=DBL) :: vnucprim(1:nprim,1:nprim)
REAL(KIND=DBL) :: vnuccont(1:ngauss,1:ngauss)
REAL(KIND=DBL) :: tempcont(1:ngauss,1:nprim)

! Temp variables
INTEGER :: i, j, jj, k, info, nbasis
INTEGER :: oswitch, rfac
REAL(KIND=DBL) :: anorm, bnorm
REAL(KIND=DBL) :: sum, coefac, exposum
REAL(KIND=DBL) :: zoff
REAL(KIND=DBL) :: partitioning
REAL(KIND=DBL) :: parfact

! Call quadrature for the Gauss-Lobotto rule
nbasis = (norder-1)*nelements-1

CALL gauss_rule(0, norder, gausspts, gausswts, work, 0.0D0, 0.0D0, 'B', info)

! Read in FEM sizes
OPEN(UNIT=1005, FILE='FEM.in', STATUS='OLD', ACTION='READ')
DO i=1, nelements
   READ(1005, *) elementsize(i)
ENDDO
CLOSE(1005)

! Make real quadrature grid
ALLOCATE(pts(1:nbasis), wts(1:nbasis))
ALLOCATE(dvrpts(1:(norder-1)*nelements+1), dvrwts(1:(norder-1)*nelements+1))

sum = 0.0d0
j = 0
jj = 0

DO i=1, nelements
   DO k=1, norder
      jj=jj+1

      ! Counter to omit the first endpt of each element except the first
      IF( (k /= 1) .OR. (i == 1)) THEN
         j = j + 1
      ENDIF

      dvrpts(j) = elementsize(i) * 0.5d0 * (gausspts(k) + 1.0d0) + sum
      dvrwts(j) = dvrwts(j) + gausswts(k) * elementsize(i) * 0.5d0

   ENDDO
   sum = sum + elementsize(i)
ENDDO


DO i=2, nbasis+1
   pts(i-1) = dvrpts(i)
   wts(i-1) = dvrwts(i)
ENDDO


! Read in exponets
! Open the file!
OPEN(UNIT=10, FILE="expos.dat", STATUS="OLD", ACTION="READ")

DO i=1, nprim
   READ(10,*) zoff, sym(i), expos(i)
ENDDO
CLOSE(10)


! Loop over gaussian basis to create prim matrix
! Check if the functions are orthogonal.  If so, the matrix element is zero
! If not, raise the exponent on r to 1+the angular r's.
DO i=1, nprim
   CALL coef(sym(i), expos(i), anorm)
   DO j=i, nprim
      CALL ortho(sym(i), sym(j), oswitch, coefac, rfac)
      sum = 0.0d0

      IF(oswitch .NE. 0) THEN
         CALL coef(sym(j), expos(j), bnorm)
         exposum = expos(i)+expos(j)
         DO k=norder, 2*norder-2
            parfact = partitioning(pts(k), pts(norder-1), pts(2*norder-2))
            sum = sum + wts(k)*EXP(-exposum*pts(k)**2)*(pts(k))**(1+rfac)*parfact
         ENDDO
         DO k=2*norder-1, nbasis
            sum = sum + wts(k)*EXP(-exposum*pts(k)**2)*(pts(k))**(1+rfac)
         ENDDO
         
      ENDIF
      vnucprim(i,j) = -sum*anorm*bnorm*PI*nucharge*coefac
      IF(ABS(vnucprim(i,j)) .LT. tol) vnucprim(i,j) = 0.0d0
      IF(j>i) vnucprim(j,i) = vnucprim(i,j)
   ENDDO
ENDDO

! xform to contracted
IF (ngauss .EQ. nprim) THEN
   xform = 0.0d0
   DO i=1, ngauss
      xform(i,i) = 1.0d0
   ENDDO
   
ELSE
   ! Call for contraction matrix
   OPEN(UNIT=95, FILE='xform.dat', STATUS='OLD', ACTION='READ')
   DO i=1, ngauss
      READ(95, *) (xform(i,j),j=1,nprim)
   ENDDO
   CLOSE(95)
ENDIF

tempcont = MATMUL(xform,vnucprim)
vnuccont = MATMUL(tempcont,TRANSPOSE(xform))

! Write out matrix
OPEN(UNIT=85, FILE='vnucprim.dat', STATUS='UNKNOWN', ACTION='WRITE')
DO i=1, nprim
   WRITE(85, *) (vnucprim(i,j), j=1, nprim)
ENDDO
CLOSE(85)

! Write out matrix
OPEN(UNIT=85, FILE='vnucpart.dat', STATUS='UNKNOWN', ACTION='WRITE')
DO i=1, ngauss
   WRITE(85, *) (vnuccont(i,j), j=1, ngauss)
ENDDO
CLOSE(85)


DEALLOCATE(pts,wts,dvrpts,dvrwts)
ENDSUBROUTINE vnucpart
