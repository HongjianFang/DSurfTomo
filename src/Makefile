CMD = SurfTomo
FC = gfortran
FFLAGS  = -O3 -ffixed-line-length-none -ffloat-store\
           -W  -fbounds-check -m64 -mcmodel=medium
F90SRCS = lsmrDataModule.f90 lsmrblasInterface.f90  \
          lsmrblas.f90 lsmrModule.f90 delsph.f90\
	  forwardstep.f90 forwardtrans.f90 split.f90 merge1.f90\
	  invtrans3d.f90 inversetrans.f90 aprod.f90\
	  haar.f90 waveletD8.f90 inversestep.f90\
	  gaussian.f90 main.f90
FSRCS =  surfdisp96.f  
OBJS = $(F90SRCS:%.f90=%.o) $(FSRCS:%.f=%.o) CalSurfG.o wavelettrans3domp.o
all:$(CMD)
$(CMD):$(OBJS)
	$(FC) -fopenmp $^ -o $@
CalSurfG.o:CalSurfG.f90
	$(FC) -fopenmp $(FFLAGS) -c $<  -o $@
wavelettrans3domp.o: wavelettrans3domp.f90
	$(FC) -fopenmp $(FFLAGS) -c $<  -o $@
%.o: %.f90
	$(FC) $(FFLAGS) -c $(@F:.o=.f90) -o $@
%.o: %.f
	$(FC) $(FFLAGS) -c $(@F:.o=.f) -o $@
clean:
	rm *.o *.mod $(CMD)
