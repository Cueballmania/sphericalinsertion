subroutine lagmatrix (nbasis, n, alfa, x, cx, lamb )
!
!*******************************************************************************
!
!! LAGUERRE_GEN evaluates the generalized Laguerre polynomials at X.
!
!
!  Differential equation:
!
!    X * Y'' + (ALFA+1-X) * Y' + N * Y = 0
!
!  Recursion:
!
!    L(0,ALFA) = 1
!    L(1,ALFA) = 1+ALFA-X
!
!    L(N,ALFA) = ( (2*N-1+ALFA-X) * L(N-1,ALFA) - (N-1+ALFA) * L(N-2,ALFA) ) / N
!
!  Restrictions:
!
!    ALFA > -1
!
!  Special values:
!
!    For ALFA = 0, the generalized Laguerre polynomial becomes the
!    Laguerre polynomial.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, real ALFA, a parameter which is part of the definition of
!    the generalized Laguerre polynomials.  It must be greater than -1.
!
!    Input, real X, the point at which the polynomials are to be
!    evaluated.
!
!    Output, real CX(0:N), the values of the first N+1 generalized
!    Laguerre polynomials at the point X.
!
  implicit none
!
  integer nbasis,n,j
!
  real(KIND=8) alfa, lamb
  COMPLEX*16 cx(1:nbasis,0:n)
  integer i
  COMPLEX*16 x(1:nbasis), lx
!
!
  if ( alfa <= -1.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LAGUERRE_GEN - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Input value of ALFA is ', alfa
    write ( *, '(a)' ) '  but ALFA must be greater than -1.'
    stop
  end if
 
  if ( n < 0 ) then
    return
  end if

!
!
WRITE(*,*) lamb
DO j=1, nbasis
  lx=lamb*x(j)
  cx(j,0) = ZEXP(-lx/2.0d0)

  if ( n == 0 ) then
    return
  end if

  cx(j,1) = ((1.0d0,0.0d0) + alfa - lx)*ZEXP(-lx/2.0d0)
!
  do i = 2, n
    cx(j,i) = ( (  (2.0d0,0.0d0) * i - (1.0d0,0.0d0)  + alfa - lx ) * cx(j,i-1) &
      - (  i - (1.0d0,0.0d0)  + alfa ) * cx(j,i-2) ) / real ( i )
  end do

  DO i=0,n
     cx(j,i) = DSQRT(lamb/((i+1.0d0)*(i+2.0d0)))*lx*cx(j,i)
  ENDDO
ENDDO
  return
end
