BINDIR = ..

FC = gfortran
FFLAGS = -Wall -std=f95    -fPIC

#LIB = -llapack -lblas
LIB = -Wl,-framework -Wl,veclib

.SUFFIXES : .o .f90 .f95 .f08 .f

OBJS =  mainscatter.f90 gaussquad.f95 DERV_DVR.f90 coulpot_insert.f90 lagpot.f90 \
		engspec.f90 gaussmat.f90 gaussmat_sch.f90 SmxInv.f90 kenginsert.f90 lagkenet.f90\
		projection.f90 SVD_orthogonalization.f95 algama.f dvr.f90 jwkb.f coul90.f \
		readmesa.f90 stos.f90 eigentest.f90 calcheck.f90 gausspot.f90 potential_insert.f90
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
