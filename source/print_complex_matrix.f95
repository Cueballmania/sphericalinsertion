SUBROUTINE print_complex_matrix (matrix,n)
!
!  generic print of a complex matrix
!  in f90 or f95 I generally allocate the matrix
!  at exactly the problem size, hence the simplicity here 
!
USE inputvariables
IMPLICIT NONE
INTEGER :: i,j,n
COMPLEX*16, DIMENSION(n,n)  :: matrix
!
! print by rows 
!
WRITE(72,'(" matrix by rows ")')
DO i=1,n
 WRITE(72,22) i,(matrix(i,j),j=1,n)
ENDDO
WRITE(72,'("  ")')
22 FORMAT('(i4,2x,5(e12.4,",",e12.4,2x),/,(6x,5(e12.4,",",e12.4,2x)))')
!
RETURN
END


