CMD = DSurfTomo
FC = gfortran
FFLAGS  = -O -ffixed-line-length-none -ffloat-store\
           -fbounds-check -m64 -mcmodel=medium
F90SRCS = lsmrDataModule.f90 lsmrblasInterface.f90\
          lsmrblas.f90 lsmrModule.f90 delsph.f90\
	  aprod.f90 gaussian.f90 voronoiproj.f90
FSRCS =  surfdisp96.f slarnv.f slaruv.f
OBJS = $(F90SRCS:%.f90=%.o) $(FSRCS:%.f=%.o) CalSurfG.o main.o
all:$(CMD)
$(CMD):$(OBJS)
	$(FC) -fopenmp $^ -o $@
CalSurfG.o:CalSurfG.f90
	$(FC) -fopenmp $(FFLAGS) -c $<  -o $@
main.o:main.f90
	$(FC) -fopenmp $(FFLAGS) -c $<  -o $@
%.o: %.f90
	$(FC) $(FFLAGS) -c $(@F:.o=.f90) -o $@
%.o: %.f
	$(FC) $(FFLAGS) -c $(@F:.o=.f) -o $@
clean:
	rm *.o *.mod $(CMD)
