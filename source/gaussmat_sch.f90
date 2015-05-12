SUBROUTINE gaussmat_sch(gaussexp, overlaps, gridpts, lagmat)
!
!Orthogonal gaussian eqvaluation on the DVR grid
!
USE inputvariables
IMPLICIT NONE
!
INTEGER :: i, j, k
COMPLEX*16 :: dvrpt, sum
REAL*8 :: gexp
REAL*8 :: gaussexp(1:numgauss),overlaps(1:numgauss, 1:numgauss)
COMPLEX*16 :: gridpts(1:nbasis), gaussdvr(1:nbasis,1:numgauss),lagmat(1:nbasis,1:numgauss)
!
DO i=1, nbasis
   dvrpt = gridpts(i)
   DO j=1, numgauss
      gexp = gaussexp(j)
      gaussdvr(i,j) = -(2.0d0*gexp/pi)**(.75)*ZEXP(-gexp*dvrpt*dvrpt)/dvrpt
   ENDDO
   WRITE(9232,*) DREAL(dvrpt), (gaussexp(j),DREAL(gaussdvr(i,j)),j=1,numgauss)
ENDDO
!
lagmat = MATMUL(gaussdvr,overlaps)
!
DO i=1, nbasis
   WRITE(2020,*) (lagmat(i,j),j=1,numgauss)
ENDDO
!
!DO i=1, nbasis
!   dvrpt = gridpts(i)
!   DO j=1, numgauss
!      sum=czero
!      DO k=1, numgauss
!         gexp = gaussexp(j)
!         sum= sum + overlaps(j,k)*(2.0d0*gexp/pi)**(.75)*ZEXP(-gexp*dvrpt*dvrpt)
!      ENDDO
!      lagmat(i,j) = sum
!   ENDDO
!ENDDO
END SUBROUTINE
