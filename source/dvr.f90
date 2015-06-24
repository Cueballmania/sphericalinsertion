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
COMPLEX(KIND=DBL) :: biggrid(1:norder*nelements)

! GaussQ variables
INTEGER :: icode = 0                 ! Legendre Quadrature
REAL(KIND=DBL) :: rpts(1:norder), lpts(1:norder)      ! Quadrature points
REAL(KIND=DBL) :: rwts(1:norder), lwts(1:norder)      ! Quadrature weights
REAL(KIND=DBL) :: work(1:norder)     ! Work array
CHARACTER(LEN=1) :: rendpts='R', lendpts='B'       ! Include both endpoints - Gauss-Lobotto
INTEGER :: info                      ! Error check

! Element size array
REAL(KIND=DBL), DIMENSION(1:nelements) :: elementsize

! Local Variables
INTEGER :: i, j, k, l, kk, jj
INTEGER :: ndvr
REAL(KIND=DBL) :: rderiv(1:norder, 1:norder)       ! First derivative matrix of the DVR functions on [-1:1] with amplitude of 1
REAL(KIND=DBL) :: lderiv(1:norder, 1:norder)
COMPLEX(KIND=DBL) :: sum, weightfac
COMPLEX(KIND=DBL) :: largederv(1:nbasis+1, 1:norder*nelements)
COMPLEX(KIND=DBL) :: smallerderv(1:nbasis,1:norder*nelements)



! GaussQ call for Radau quadrature points on [-1:1]
CALL gauss_rule(icode, norder, rpts, rwts, work, 0.0D0, 0.0D0, rendpts, info)
IF(info .NE. 0) THEN
     WRITE(*,*) 'gauss_rule error', info
     STOP
ENDIF


! GaussQ call for Lobotto quadrature points on [-1:1]
CALL gauss_rule(icode, norder, lpts, lwts, work, 0.0D0, 0.0D0, lendpts, info)
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

! number of points with the last endpoint
ndvr = (norder-1)*nelements + 1


! Assemble the grid
! First part of the grid is Gauss-Radau quadrature points
sum = (0.0d0,0.0d0)
j = 0
jj = 0
i=1

weights = (0.0d0,0.0d0)
weightfac = (1.0d0,0.0d0)

! Check to see if we're on the compelx contour
IF( i > (nelements - celements) ) THEN
   weightfac = EXP((0.0d0,1.0d0) * (angle*PI/180.0d0))
ENDIF

DO k=1, norder
   jj = jj + 1
   j = j + 1

   dvrpts(j) = weightfac * elementsize(i) * 0.5d0 * (rpts(k) + 1.0d0) + sum
   dvrwts(j) = dvrwts(j) + rwts(k) * elementsize(i) * 0.5d0 * weightfac*dvrpts(j)*dvrpts(j)
   biggrid(jj) = dvrpts(j)
   bigweights(jj) = rwts(k) * elementsize(i) * 0.5d0 * weightfac

ENDDO

sum = sum + elementsize(i) * weightfac






DO i=2, nelements
   weightfac = (1.0d0,0.0d0)
  
   ! Check to see if we're on the compelx contour
   IF( i > (nelements - celements) ) THEN
      weightfac = EXP((0.0d0,1.0d0) * (angle*PI/180.0d0))
   ENDIF

   DO k=1, norder
      jj=jj+1

      ! Counter to omit the first endpt of each element except the first
      IF( (k /= 1) ) THEN
         j = j + 1
      ENDIF

      dvrpts(j) = weightfac * elementsize(i) * 0.5d0 * (lpts(k) + 1.0d0) + sum
      dvrwts(j) = dvrwts(j) + lwts(k) * elementsize(i) * 0.5d0 * weightfac*dvrpts(j)*dvrpts(j)
      biggrid(jj) = dvrpts(j)
      bigweights(jj) = lwts(k) * elementsize(i) * 0.5d0 * weightfac

   ENDDO

   sum = sum + elementsize(i) * weightfac
ENDDO


! Shrink the grid and the weights by dropping the last point
DO i=1, nbasis
   grid(i) = dvrpts(i)
   weights(i) = dvrwts(i)
ENDDO

OPEN (unit=300, file="dvrgrid.out", status='unknown')
DO i=1, nbasis
   WRITE(300,'(1x,8000ES17.5)') DREAL(grid(i)), DIMAG(grid(i)), weights(i)
ENDDO
CLOSE(300)






! Create Radau first derivative matrix for the quadrature - not normalized or scaled to the actual grid
DO i=1 , norder
   DO j=1, norder
      IF(i /= j) THEN
         rderiv(i,j) = 1.0d0/(rpts(i)-rpts(j))
         DO k=1, norder
            IF( k /= i .AND. k /= j) THEN
               rderiv(i,j) = rderiv(i,j)*(rpts(j)-rpts(k))/(rpts(i)-rpts(k))
            ENDIF
         ENDDO
      ELSE
         rderiv(i,i)=0.0d0
         DO k=1, norder
            IF(k /= i) THEN
               rderiv(i,i) = rderiv(i,i) + 1.0d0/(rpts(i)-rpts(k))
            ENDIF
         ENDDO
      ENDIF
   ENDDO
ENDDO

! Create Lobotto first derivative matrix for the quadrature - not normalized or scaled to the actual grid
DO i=1 , norder
   DO j=1, norder
      IF(i /= j) THEN
         lderiv(i,j) = 1.0d0/(lpts(i)-lpts(j))
         DO k=1, norder
            IF( k /= i .AND. k /= j) THEN
               lderiv(i,j) = lderiv(i,j)*(lpts(j)-lpts(k))/(lpts(i)-lpts(k))
            ENDIF
         ENDDO
      ELSE
         lderiv(i,i)=0.0d0
         DO k=1, norder
            IF(k /= i) THEN
               lderiv(i,i) = lderiv(i,i) + 1.0d0/(lpts(i)-lpts(k))
            ENDIF
         ENDDO
      ENDIF
   ENDDO
ENDDO






! Assemble the KE matrix
! Create the matrix of each dvr function at EVERY dvrpoint
largederv = (0.0d0,0.0d0)

i=1

weightfac = (1.0d0,0.0d0)
  
! Check to see if we're on the compelx contour
IF( i > (nelements - celements) ) THEN
   weightfac = EXP((0.0d0,1.0d0) * (angle*PI/180.0d0))
ENDIF

DO j=1, norder
   l = (i-1)*(norder-1)+j

   DO k=1, norder
      kk = (i-1)*norder+k    

      largederv(l,kk) = rderiv(j,k)*2.0d0/elementsize(i)/weightfac
   ENDDO
ENDDO




DO i=2, nelements
   weightfac = (1.0d0,0.0d0)
  
   ! Check to see if we're on the compelx contour
   IF( i > (nelements - celements) ) THEN
      weightfac = EXP((0.0d0,1.0d0) * (angle*PI/180.0d0))
   ENDIF

   DO j=1, norder
      l = (i-1)*(norder-1)+j

      DO k=1, norder
         kk = (i-1)*norder+k    

         largederv(l,kk) = lderiv(j,k)*2.0d0/elementsize(i)/weightfac
      ENDDO
   ENDDO
ENDDO

! Normalize largederv
DO i=1, nbasis+1
   largederv(i,:) = largederv(i,:)/SQRT(dvrwts(i))
ENDDO


! Shrink to the basis size
smallerderv = (0.0d0,0.0d0)
DO i=1, nbasis
   smallerderv(i,:) = largederv(i,:)
ENDDO

! Create the KE matrix
DO i=1, nbasis
   DO j=1, nbasis
      sum=(0.0d0,0.0d0)
      DO k=1, norder*nelements
         sum = sum + smallerderv(i,k)*smallerderv(j,k)*bigweights(k)*biggrid(k)*biggrid(k)
      ENDDO
      kinetic(i,j) = sum*0.5d0
   ENDDO
ENDDO



END SUBROUTINE
