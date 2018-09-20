FC      = $(shell printenv FC)
ifeq ($(FC),)
	FC = gfortran
endif
CFLAGS  = -g -O3 -march=core-avx2 -fopenmp
CFLAGS += -D_DEBUG

all:
	$(FC) $(CFLAGS) 1d_1st.F90

clean:
	rm -f a.out *.mod result.dat
