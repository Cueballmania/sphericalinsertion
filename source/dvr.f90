SUBROUTINE DVR(nelements, celements, norder, nbasis, angle, grid, weights, kinetic)
!SUBROUTINE DVR(nelements, celements, norder, nbasis, angle)
! This subroutine creature a FEDVR grid with either no scaling, straight complex
! scaling or exterior complex scaling
! It returns the grid, weights and the kinetic energy matrix only in r -- no centrifugal term
!
! This subroutine reades elements.inp for starting and endpoints of each element
!
USE gaussquad
IMPLICIT NONE

! Declare double precision
INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=13,r=200)
REAL(KIND=DBL), PARAMETER :: PI = 3.141592653589

! Input variables
INTEGER, INTENT(IN) :: nelements   ! number of total elements
INTEGER, INTENT(IN) :: celements   ! number of complex elements at the tail of the grid
INTEGER, INTENT(IN) :: norder      ! order of the quadrature in each element
INTEGER, INTENT(IN) :: nbasis      ! number of gridpts in the basis

REAL(KIND=DBL), INTENT(IN) :: angle     ! angle of the complex scaling

! Output variables
COMPLEX(KIND=DBL), INTENT(OUT) :: grid(1:nbasis)                 ! grid of points
COMPLEX(KIND=DBL), INTENT(OUT) :: weights(1:nbasis)              ! DVR weights
COMPLEX(KIND=DBL), INTENT(OUT) :: kinetic(1:nbasis,1:nbasis)     ! KE matrix

! Internal Variables - basis + endpoints
COMPLEX(KIND=DBL) :: dvrpts(1:(norder-1)*nelements+1) ! DVR points
COMPLEX(KIND=DBL) :: dvrwts(1:(norder-1)*nelements+1) ! DVR weights
COMPLEX(KIND=DBL) :: bigweights(1:norder*nelements)

! GaussQ variables
INTEGER :: icode = 0                 ! Legendre Quadrature
REAL(KIND=DBL) :: pts(1:norder)      ! Quadrature points
REAL(KIND=DBL) :: wts(1:norder)      ! Quadrature weights
REAL(KIND=DBL) :: work(1:norder)     ! Work array
CHARACTER(LEN=1) :: endpts='B'       ! Include both endpoints - Gauss-Lobotto
INTEGER :: info                      ! Error check

! Element size array
REAL(KIND=DBL), DIMENSION(1:nelements) :: elementsize

! Local Variables
INTEGER :: i, j, k, l, kk, jj
INTEGER :: ndvr
REAL(KIND=DBL) :: deriv(1:norder, 1:norder)       ! First derivative matrix of the DVR functions on [-1:1] with amplitude of 1
COMPLEX(KIND=DBL) :: sum, weightfac
COMPLEX(KIND=DBL) :: largederv(1:nbasis+2, 1:norder*nelements)
COMPLEX(KIND=DBL) :: smallerderv(1:nbasis,1:norder*nelements)

! GaussQ call for quadrature points on [-1:1]
CALL gauss_rule(icode, norder, pts, wts, work, 0.0D0, 0.0D0, endpts, info)
IF(info .NE. 0) THEN
     WRITE(*,*) 'gauss_rule error', info
     STOP
ENDIF


! Read in FEM sizes
OPEN(UNIT=1005, FILE='FEM.in', STATUS='OLD', ACTION='READ')
DO i=1, nelements
   READ(1005, *) elementsize(i)
ENDDO
CLOSE(1005)

ndvr = (norder-1)*nelements + 1

! Assemble the grid
sum = (0.0d0,0.0d0)
j = 0
jj = 0

weights = (0.0d0,0.0d0)

DO i=1, nelements
   weightfac = (1.0d0,0.0d0)
  
   ! Check to see if we're on the compelx contour
   IF( i > (nelements - celements) ) THEN
      weightfac = EXP((0.0d0,1.0d0) * (angle*PI/180.0d0))
   ENDIF

   DO k=1, norder
      jj=jj+1

      ! Counter to omit the first endpt of each element except the first
      IF( (k /= 1) .OR. (i == 1)) THEN
         j = j + 1
      ENDIF

      dvrpts(j) = weightfac * elementsize(i) * 0.5d0 * (pts(k) + 1.0d0) + sum
      dvrwts(j) = dvrwts(j) + wts(k) * elementsize(i) * 0.5d0 * weightfac
      bigweights(jj) = wts(k) * elementsize(i) * 0.5d0 * weightfac

   ENDDO

   sum = sum + elementsize(i) * weightfac
ENDDO


! Shrink the grid and the weights to only have the basis functions
DO i=2, nbasis+1
   grid(i-1) = dvrpts(i)
   weights(i-1) = dvrwts(i)
ENDDO

OPEN (unit=300, file="dvrgrid.out", status='unknown')
DO i=1, nbasis
   WRITE(300,'(1x,8000ES17.5)') DREAL(grid(i)), DIMAG(grid(i)), weights(i)
ENDDO
CLOSE(300)

! Create first derivative matrix for the quadrature - not normalized or scaled to the actual grid
DO i=1 , norder
   DO j=1, norder
      IF(i /= j) THEN
         deriv(i,j) = 1.0d0/(pts(i)-pts(j))
         DO k=1, norder
            IF( k /= i .AND. k /= j) THEN
               deriv(i,j) = deriv(i,j)*(pts(j)-pts(k))/(pts(i)-pts(k))
            ENDIF
         ENDDO
      ELSE
         deriv(i,i)=0.0d0
         DO k=1, norder
            IF(k /= i) THEN
               deriv(i,i) = deriv(i,i) + 1.0d0/(pts(i)-pts(k))
            ENDIF
         ENDDO
      ENDIF
   ENDDO
ENDDO

! Assemble the KE matrix
! Create the matrix of each dvr function at EVERY dvrpoint
largederv = (0.0d0,0.0d0)

DO i=1, nelements
   weightfac = (1.0d0,0.0d0)
  
   ! Check to see if we're on the compelx contour
   IF( i > (nelements - celements) ) THEN
      weightfac = EXP((0.0d0,1.0d0) * (angle*PI/180.0d0))
   ENDIF

   DO j=1, norder
      l = (i-1)*(norder-1)+j

      DO k=1, norder
         kk = (i-1)*norder+k    

         largederv(l,kk) = deriv(j,k)*2.0d0/elementsize(i)/weightfac
      ENDDO
   ENDDO
ENDDO

! Normalize largederv
DO i=1, nbasis+2
   largederv(i,:) = largederv(i,:)/SQRT(dvrwts(i))
ENDDO


! Shrink to the basis size
smallerderv = (0.0d0,0.0d0)
DO i=2, nbasis+1
   smallerderv(i-1,:) = largederv(i,:)
ENDDO

! Create the KE matrix
DO i=1, nbasis
   DO j=1, nbasis
      sum=(0.0d0,0.0d0)
      DO k=1, norder*nelements
         sum = sum + smallerderv(i,k)*smallerderv(j,k)*bigweights(k)
      ENDDO
      kinetic(i,j) = sum*0.5d0
   ENDDO
ENDDO


!!$OPEN (unit=9996, file="dvrKE.out", status='unknown')
!!$DO i=1, nbasis
!!$   WRITE(9996,'(1x,8000ES17.5)') (DREAL(kinetic(i,j)),DIMAG(kinetic(i,j)) ,j=1,nbasis)
!!$ENDDO
!!$CLOSE(9996)

END SUBROUTINE
