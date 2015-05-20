BINDIR = ..

FC = gfortran
FFLAGS = -Wall -std=f95    -fPIC

#LIB = -llapack -lblas
LIB = -Wl,-framework -Wl,veclib

.SUFFIXES : .o .f90 .f95 .f08 .f

OBJS =  mainscatter.f90 gaussquad.f95 engspec.f90 gaussmat.f90 algama.f \
		SVD_orthogonalization.f95 dvr.f90 readmesa.f90 stos.f90 eigentest.f90 \
		calcheck.f90 gausspot.f90 insertion.f90 partitioning.f90 vnucpart.f90\
                ortho.f90 coef.f90
.f.0:
	$(FC) -c $(FFLAGS) $<

.f08.o:
	$(FC) -c $(FFLAGS) $<

.f90.o:
	$(FC) -c $(FFLAGS) $<

.f95.o:
	$(FC) -c $(FFLAGS) $<

xscattering: $(OBJS) gaussquad.mod
	$(FC) $(OBJS) -o $(BINDIR)/$@ $(LIB)
