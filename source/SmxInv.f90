SUBROUTINE SmxInv (a,n)
!USE Numeric_Kinds_Module
!
!
! >>> SmxInv    -- computes in-place inverse of symmetric matrix
!
! *** Calling Sequence:   CALL SmxInv ( a,n )
!
! === Latest Revision:  |  |  set upper half too
!                       | 921104 at 21:14:57 |  lifted from GAP routine
!     revised by rbw
! --- on entry the calling routine supplies --
!
!     n         -- size of matrix to be inverted
!     a(i,j)    -- matrix to be inverted
!
! --- on exit the routine returns --
!
!     a(i,j)    -- the inverted matrix
!     nfp       -- the number of floating point operations performed
!
! --- information from the following common blocks is used --
!
!
!
!----------------------------- local declarations
!
IMPLICIT NONE
INTEGER n, nfp, i, im1, j, np1, k, kp1, ii, ip1, ifp
REAL(KIND=8) a(n,n), den, fac, c
!
!---------------------------- execution begins here --------------------
!
IF(n.eq.1) THEN
    a(1,1) = 1d0/a(1,1)
    RETURN
END IF
ifp = 0
!
!       --- this code from GAP routine of same name ---
!  -- construct upper half of symmetric matrix from the lower---


do i = 2, n
   im1 = i-1
   do j = 1, im1
      a(i,j) = a(j,i)
   ENDDO
ENDDO
!
!  -- begin inversion
!
    np1 = n + 1
    do 40 k = 1, n
    kp1  = k + 1
    den  = 1d0/a(k,k)
    a(k,k) = 1d0
ifp = ifp + 1

    do 20 j = 1,n
        a(j,k) = a(j,k)*den
20      CONTINUE
    ifp = ifp + n
    IF(k.eq.n) GOTO 40

    do 30 i = kp1, n
        fac = a(k,i)
        a(k,i) = 0d0
        do 31 j = 1,n
        a(j,i) = a(j,i) - fac*a(j,k)
31        CONTINUE
        ifp = ifp + n
30      CONTINUE

40    CONTINUE
!
    do 60 ii = 2, n
    i = np1 - ii
    ip1 = i+1
    do 61 j = 1, i
        c = a(j,i)
        do 50 k = ip1, n
        c = c - a(k,i)*a(j,k)
50        CONTINUE
        a(j,i) = c
    ifp = ifp + n - i
61      CONTINUE
60    CONTINUE
!
!  -- construct lower half of inverted symmetrix matrix from upper--
!
do j = 2, n
   do i = 1, j-1
      a(j,i) = a(i,j)
   ENDDO
ENDDO
nfp = ifp
RETURN
END SUBROUTINE SmxInv
