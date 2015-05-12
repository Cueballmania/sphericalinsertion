!This subroutine generates the kinetic energy matrix
!It is block diagonal and uses calls to the subroutine
!basisprime which generates the derivatives of the basis
!functions
SUBROUTINE keng(kinetic, mesh, weights,gridweights)
USE inputvariables
IMPLICIT NONE
!
INTEGER :: i, j, k, l
COMPLEX*16 :: mesh(0:numelements-1,0:norder-1)
COMPLEX*16 :: kinetic(1:nbasis,1:nbasis), blockken(0:norder-1,0:norder-1)
COMPLEX*16 :: weights(0:numelements-1,0:norder-1)
COMPLEX*16 :: gridweights(1:nbasis)
COMPLEX*16 :: toolarge(0:numelements*(norder-1),0:numelements*(norder-1))
COMPLEX*16 :: integ
!
kinetic(:,:)=czero
toolarge(:,:)=czero
!
!Generates the kinetic energy matrix in each element
!
DO i=0, numelements-1
   CALL DERV_DVR(mesh(i,:),weights(i,:),blockken)
   IF(i.EQ.0) THEN
      DO l=0, norder-1
!      WRITE(773,*) (blockken(l,k), k=0, norder-1)
      ENDDO
   ENDIF

!
   DO j=0, norder-1
      DO k=0, norder-1
         integ=czero
         DO l=0,norder-1
            integ = integ + blockken(j,l)*blockken(k,l)*weights(i,l)
         ENDDO
         toolarge(i*(norder-1)+j,i*(norder-1)+k) = toolarge(i*(norder-1)+j,i*(norder-1)+k)+integ
      ENDDO
   ENDDO
ENDDO
WRITE(*,*) "Made the kinetic energy matrix too large"
!
!Make the kinetic energy matrix
!
DO i=1, nbasis
   DO j=1, nbasis
      kinetic(i,j)=toolarge(i,j)
   ENDDO
ENDDO
!
!Normalize
!
DO i=1, nbasis
   DO j=1, nbasis
      kinetic(i,j) = kinetic(i,j)/ZSQRT(gridweights(i)*gridweights(j))
   ENDDO
ENDDO
WRITE(*,*) "Made the Kinetic Energy matrix"
!
!
! Scaling the kinentic energy matrix from the derivative term
kinetic(:,:)=5D-1*kinetic(:,:)/mass
!
OPEN (unit=9996, file="dvrkinetic.out", status='unknown')
DO i=1, nbasis
   WRITE(9996,'(1x,8000ES17.5)') (DREAL(kinetic(i,j)),DIMAG(kinetic(i,j)) ,j=1,nbasis)
ENDDO
CLOSE(9996)
!
END SUBROUTINE
