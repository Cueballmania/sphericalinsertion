!This Subroutine will create an exterior complex scaling mesh
!given the number of elements, the order of Lobatto quadrature,
!the point where complex scaling begins and the angle it takes.
!There will be a seperate input to specify if one wants more than
!one element on the complex contour.
SUBROUTINE complexgrid(mesh, weights)
USE gaussquad
USE inputvariables
IMPLICIT NONE
!
REAL(KIND=8) :: startpt, endpt, elementsize
INTEGER :: icode=0, info, i, j
REAL(KIND=8) :: pts(0:norder-1), tempw(0:norder-1), work(0:norder-1)
COMPLEX*16 :: mesh(0:numelements-1,0:norder-1)
COMPLEX*16 :: weights(0:numelements-1,0:norder-1)
!
OPEN(unit=1002, file='ptsweights.out', status='unknown')
OPEN(unit=1003, file='fullgrid.out', status='unknown')
OPEN(unit=1005, file='elements.inp', status='unknown')
!
!Calls the gaussquad.f95 to set up the quadrature points.
!icode=0 for Gauss-Legendre quadrature
!norder is the order of the quandrature
!pts is the array to hold all of the quadrature points
!weights is the array to hold the weights of the quadrature points
!start and stop are dummy variables to ignore
!'B' specifies quadrature points including the Both endpoints
!info is the errorvariable
!
CALL gauss_rule(icode, norder, pts, tempw, work, 0.0D0, 0.0D0, 'B', info)
!
IF(info .NE. 0) THEN
     WRITE(*,*) 'gauss_rule error', info
     STOP
ENDIF
!
!Sets up the quadrature points with the exterior complex scaling.
!If customgrid, sets spacing approriately.
!If not, divides space into evenly partitianed elements
!
READ(1005,*) startpt, endpt
elementsize=endpt-startpt
!
!GENERATES the real quadrature points
!
DO i=0, relements-1
    DO j=0, norder-1
        weights(i,j)=0.5D0*elementsize*tempw(j)
        mesh(i,j)=5.0D-1*(pts(j)*elementsize+startpt+endpt)
        IF(DREAL(mesh(i,j)) < TOLERANCE) mesh(i,j) = czero
    ENDDO
    IF(i.NE.numelements-1) THEN
       READ(1005,*) startpt, endpt
       elementsize=endpt-startpt
    ENDIF
ENDDO
!
!GENERATES complex quadrature points
!
contourpt=startpt
!
DO i=0,celements-1
    DO j=0, norder-1
        weights(i+relements,j)=0.5d0*elementsize*tempw(j)*contourphase
        mesh(i+relements,j)=contourpt+(5.0D-1*(pts(j)*elementsize+startpt+endpt)-contourpt)*contourphase
    ENDDO
    IF(i.NE.celements-1) THEN
       READ(1005,*) startpt, endpt
       elementsize=endpt-startpt
    ENDIF
ENDDO
gridend = endpt
CLOSE(1005)
!
!Writes the quadrature points and weights out
!
DO i=0, norder-1
   WRITE(1002,*) pts(i), tempw(i)
ENDDO
!
CLOSE(1002)
!
!Writes out the exterior complex scaling grid
!
DO i=0,numelements-1
   DO j=0, norder-1
   WRITE(1003,*)  DREAL(mesh(i,j)), DIMAG(mesh(i,j)), DREAL(weights(i,j)), DIMAG(weights(i,j))
   ENDDO
ENDDO
!
ClOSE(1003)
END SUBROUTINE
