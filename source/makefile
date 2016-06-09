BINDIR = ..

FC = gfortran
FFLAGS = -lstdc++

#LIB = -llapack -lblas
LIB = -Wl,-framework -Wl,veclib

.SUFFIXES : .o .f90 .f95 .f08 .f

OBJS =  mainscatter.f90 gaussquad.f95 engspec.f90 gaussmat.f90 algama.f \
		SVD_orthogonalization.f95 dvr.f90 readmesa.f90 stos.f90 eigentest.f90 \
		calcheck.f90 gausspot.f90 insertion.f90 bumppart.f90 beckepart.f90 \
                fortran_faddeeva.o faddeeva.o direct_coul.f90
.f.0:
	$(FC) -c $(FFLAGS) $<

.f08.o:
	$(FC) -c $(FFLAGS) $<

.f90.o:
	$(FC) -c $(FFLAGS) $<

.f95.o:
	$(FC) -c $(FFLAGS) $<

xscattering: $(OBJS) gaussquad.mod
	 $(FC) $(FFLAGS) $(OBJS) -o $(BINDIR)/$@ $(LIB)
