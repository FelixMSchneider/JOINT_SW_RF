FC = gfortran
FFLAGS  = -O -ffixed-line-length-none -ffloat-store\
           -fbounds-check -m64 -mcmodel=medium

all:
	$(FC) $(FFLAGS) -c  surfdisp96.f
	gfortran -c igetmod.f
	gfortran -o call_surf call_surf.f surfdisp96.o igetmod.o
	rm *.o
run:
	python get_Disp.py
