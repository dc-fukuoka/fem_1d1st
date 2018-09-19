FC = $(shell printenv FC)
ifeq ($(FC),)
	FC = gfortran
endif
all:
	$(FC) -D_DEBUG -g -O3 -march=core-avx2 -fopenmp 1d_1st.F90

clean:
	rm -f a.out *.mod result.dat
