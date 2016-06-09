! This subroutine is the hard-coded way to calculate the direct Coulomb operator for Beryllium using 5-s functions
! that make up the two occupied orbitals.
! The first two s-functions is made up of a contraction of 6 and 2 functions respectively.

SUBROUTINE direct_coul(grid,pot)
USE inputvariables
IMPLICIT NONE

! Input of the grid and the resulting direct operator
COMPLEX(KIND=DBL), INTENT(IN) :: grid(1:nbasis)                 ! grid of points
COMPLEX(KIND=DBL), INTENT(OUT) :: pot(1:nbasis)              ! DVR weights

! The coefficents of the 5-contracted gaussians for the description of the two occupied orbitals
REAL(KIND=DBL), DIMENSION(1:2,1:5) :: c

! The transformation matrix for the 10-s fucntions to 5-s contracted functions
REAL(KIND=DBL), DIMENSION(1:5,1:10) :: a

! The Gaussian exponents for the 10 primitives
REAL(KIND=DBL), DIMENSION(1:10) :: alpha = (/3630.0, 532.3, 117.8, 32.66, 10.48, 3.668, 1.354, 0.389, 0.1502, 0.05241/)

! Indices
INTEGER :: i,j,k,l,m,n

! result of the complex error function
COMPLEX(KIND=DBL) :: fad


c = TRANSPOSE(RESHAPE((/6.26743227E-01,4.16593765E-01,3.70280025E-02,-7.87718703E-03,1.82218542E-03,1.21177112E-01,&
1.88531880E-01,-4.50329330E-02,-5.93418593E-01,-4.70972509E-01/), (/ 5,2 /)))

a = TRANSPOSE(RESHAPE((/0.0008390,0.0067350,0.0357260,0.1386350,0.3853990,0.5476880,0.0,0.0,0.0,0.0,&
0.0,0.0,0.0,0.0,0.0,0.2134060,0.8146920,0.0,0.0,0.0,0.0,0.0,&
0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.,0.0,0.0,0.0,0.0,0.0,0.0,&
0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0/), (/ 10,5 /)))

pot(:) = czero
DO n=1, nbasis
   DO i=1, 2
      DO j=1,5
         DO k=1,5
            DO l=1, 10
               DO m=1, 10
                  CALL faddeeva_erf(fad,SQRT(alpha(l)+alpha(m))*grid(n))
                  pot(n) = pot(n)+c(i,j)*c(i,k)*a(j,l)*a(k,m)*(2.0/(alpha(l)&
+alpha(m)))**(1.5)*(alpha(l)*alpha(m))**(0.75)*fad/grid(n)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

OPEN(UNIT=23, FILE='directpot.out', STATUS='UNKNOWN', ACTION='WRITE')
DO i=1, nbasis
   WRITE(23, *) REAL(grid(i)), IMAG(grid(i)), REAL(pot(i)), IMAG(pot(i))
ENDDO
CLOSE(23)

END SUBROUTINE
