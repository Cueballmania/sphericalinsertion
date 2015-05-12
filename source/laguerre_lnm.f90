subroutine laguerre_lnm ( n, m, x, cx )
!
!*******************************************************************************
!
!! LAGUERRE_LNM evaluates the associated Laguerre polynomials Lnm at X.
!
!
!  Differential equation:
!
!    X Y'' + (M+1-X) Y' + (N-M) Y = 0
!
!  First terms:
!
!    M = 0
!
!    L(0,0) =   1
!    L(1,0) =  -X    +  1
!    L(2,0) =   X**2 -  4 X     +  2
!    L(3,0) =  -X**3 +  9 X**2 -  18 X    +    6
!    L(4,0) =   X**4 - 16 X**3 +  72 X**2 -   96 X +      24
!    L(5,0) =  -X**5 + 25 X**4 - 200 X**3 +  600 X**2 -  600 x    +  120
!    L(6,0) =   X**6 - 36 X**5 + 450 X**4 - 2400 X**3 + 5400 X**2 - 4320 X + 720
!
!    M = 1
!
!    L(0,1) =    0
!    L(1,1) =   -1,
!    L(2,1) =    2 X-4,
!    L(3,1) =   -3 X**2 + 18 X - 18,
!    L(4,1) =    4 X**3 - 48 X**2 + 144 X - 96
!
!    M = 2
!
!    L(0,2) =    0
!    L(1,2) =    0,
!    L(2,2) =    2,
!    L(3,2) =   -6 X + 18,
!    L(4,2) =   12 X**2 - 96 X + 144
!
!    M = 3
!
!    L(0,3) =    0
!    L(1,3) =    0,
!    L(2,3) =    0,
!    L(3,3) =   -6,
!    L(4,3) =   24 X - 96
!
!    M = 4
!
!    L(0,4) =    0
!    L(1,4) =    0
!    L(2,4) =    0
!    L(3,4) =    0
!    L(4,4) =   24
!
!  Recursion:
!
!    if N < M:
!      L(N,M)   = 0
!    if N = M:
!      L(N,M)   = (-1)**M * M! 
!    if N = M+1:
!      L(N,M) = (M+1)*(M+1-X)*L(M,M)
!    if N = M+2 or greater:
!      L(N,M)  = -( (X+M-2*N+1)*L(N-1,M) + (N-1)*(N-1)*L(N-2,M) )*N / (N-M)
!
!  Special values:
!
!    For M = 0, the associated Laguerre polynomials are equal to the
!    Laguerre polynomials.
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
!    Input, integer M, the parameter.  M must be nonnegative.
!
!    Input, real X, the point at which the polynomials are to be evaluated.
!
!    Output, real CX(0:N), the values of the first N+1 Laguerre
!    polynomials at the point X.
!
  implicit none
!
  integer n
!
  COMPLEX*16 :: cx(0:n)
  integer i
  integer ifact
  integer m
  COMPLEX*16 :: x
!
  if ( m < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LAGUERRE_LNM - Fatal error!'
    write ( *, '(a,i6)' ) '  Input value of M = ', m
    write ( *, '(a)' ) '  but M must be nonnegative.'
    stop
  end if
 
  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LAGUERRE_LNM - Fatal error!'
    write ( *, '(a,i6)' ) '  Input value of N = ', n
    write ( *, '(a)' ) '  but N must be positive.'
    stop
  end if

  if ( n < 0 ) then
    return
  end if

  cx(:) = (0.0d0, 0.0d0)

  ifact = 1
  do i = 1, m
    ifact = - ifact * i
  end do
 
  cx(m) = DCMPLX( ifact )
  cx(m+1) = DCMPLX ( m + 1 ) * ( DCMPLX ( m + 1 ) - x ) * cx(m)

  do i = m+2, n
    cx(i) = - ( DCMPLX ( i ) * &
      ( x + DCMPLX ( m - 2 * i + 1 ) ) * cx(i-1) &
      + DCMPLX ( ( i - 1 )**2 ) * cx(i-2) ) / DCMPLX ( i - m )
  end do

  return
end
