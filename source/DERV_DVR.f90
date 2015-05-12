! This subroutine evaluates the derivate of a DVR function in an element
SUBROUTINE DERV_DVR(pts,weights,blockken)
USE inputvariables
IMPLICIT NONE
!
INTEGER :: i,j,l,k
COMPLEX*16 :: pts(0:norder-1)
COMPLEX*16 :: weights(0:norder-1)
COMPLEX*16 :: blockken(0:norder-1,0:norder-1)
!
blockken(:,:)=czero
!
DO j=0, norder-1
   DO k=0, norder-1
      IF (k .NE. j) THEN
         blockken(k,j) = cone/(pts(k)-pts(j))
         DO l=0, norder-1
            IF (l .NE. k .AND. l .NE. j) THEN
               blockken(k,j) = blockken(k,j) * (pts(j)-pts(l))/(pts(k)-pts(l))
            ENDIF
         ENDDO
      ELSE
         IF (k .EQ. 0) THEN
            blockken(k,k) = -0.5d0/weights(k)
         ELSEIF(k .EQ. norder-1) THEN
            blockken(k,k) = 0.5d0/weights(k)
         ELSE
            blockken(k,k) = czero
         ENDIF
      ENDIF
   ENDDO
ENDDO
!
ENDSUBROUTINE
