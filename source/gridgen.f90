!This subroutine generates the diagonal values of the interpolating polynomials.
!It's also the diagonal values of the basis functions noting that
!the basis functions only have non-zero values at one value.
SUBROUTINE gridgen(gridpts, mesh)
USE inputvariables
IMPLICIT NONE
!
INTEGER :: i, j, felement, nindex
COMPLEX*16 :: mesh(0:numelements-1,0:norder-1)
COMPLEX*16 :: gridpts(1:nbasis)
!
!
!Generates the basis functions
!
j=1
DO i=1, norder*numelements-1
   IF(MOD(i,norder).NE.norder-1) THEN
      nindex = MOD(i,norder)
      felement = i/norder
      gridpts(j) = mesh(felement,nindex)
      j=j+1
   ENDIF
ENDDO
!
!Output the basis
!
OPEN(unit=9994, file='grid.out', status='unknown')
!
DO i=1,nbasis
   WRITE(9994,*) DREAL(gridpts(i)), DIMAG(gridpts(i))
ENDDO
!
CLOSE(9994)
!
END SUBROUTINE
