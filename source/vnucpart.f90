SUBROUTINE vnucpart
! This subroutine utilizes the underlying quadrature
! to re-calculate the nuclear attraction
! matrix elements using the partitioning function.
! This subroutine returns nothing, but generates the file
! vnucpart.dat
!
! This routine also reads FEM.in and xform.dat

! Input variables
INTEGER, INTENT(IN) :: norder
INTEGER, INTENT(IN) :: nelements
INTEGER, INTENT(IN) :: ngauss
INTEGER, INTENT(IN) :: nprim

! Quadrature arrays
REAL(KIND=DBL) :: pts
REAL(KIND=DBL) :: wts

! Gauss basis info
REAL(KIND=DBL) :: expos(1:nprim)
CHAR(LEN=3) :: sym(1:nprim)
REAL(KIND=DBL) :: xform(1:ngauss,1:nprim)

! Matrices
REAL(KIND=DBL) :: vnucprim(1:nprim,1:nprim)
REAL(KIND=DBL) :: vnuccont(1:ngauss,1:ngauss)

! Temp variables
INTEGER :: i, j
REAL(KIND=DBL) :: sum
REAL(KIND=DBL) :: partitionin

! Call quadrature



! Make real quadrature grid



! Read in exponets


! Loop over gaussian basis to create prim matrix



! xform to contracted


! Write out matrix
