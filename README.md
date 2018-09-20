fem_1d1st
======
one dimensional first order element finite element method(Galerkin method) using BiCGSTAB method.

requirements for python:
------
* numpy
* scipy
* matplotlib
* f90nml

how to test:
------

~~~
$ vi fort.11 # adjust length, this is the right boundary, left boundary is zero
$ vi fort.7  # add values between zero and length
$ make
$ ./a.out
 e%nelements:           7
 e%nnodes:              8
 points:
   0.00000E+00   1.00000E+00   2.00000E+00   2.10000E+00   2.30000E+00   2.70000E+00   2.90000E+00   3.00000E+00
 lhs A:
   1.00000E+00   0.00000E+00   0.00000E+00   0.00000E+00   0.00000E+00   0.00000E+00   0.00000E+00   0.00000E+00
  -1.00000E+00   2.00000E+00  -1.00000E+00   0.00000E+00   0.00000E+00   0.00000E+00   0.00000E+00   0.00000E+00
   0.00000E+00  -1.00000E+00   1.10000E+01  -1.00000E+01   0.00000E+00   0.00000E+00   0.00000E+00   0.00000E+00
   0.00000E+00   0.00000E+00  -1.00000E+01   1.50000E+01  -5.00000E+00   0.00000E+00   0.00000E+00   0.00000E+00
   0.00000E+00   0.00000E+00   0.00000E+00  -5.00000E+00   7.50000E+00  -2.50000E+00   0.00000E+00   0.00000E+00
   0.00000E+00   0.00000E+00   0.00000E+00   0.00000E+00  -2.50000E+00   7.50000E+00  -5.00000E+00   0.00000E+00
   0.00000E+00   0.00000E+00   0.00000E+00   0.00000E+00   0.00000E+00  -5.00000E+00   1.50000E+01  -1.00000E+01
   0.00000E+00   0.00000E+00   0.00000E+00   0.00000E+00   0.00000E+00   0.00000E+00  -1.00000E+01   1.00000E+01

 rhs b:
   0.00000E+00
  -1.00000E+00
  -5.50000E-01
  -1.50000E-01
  -3.00000E-01
  -3.00000E-01
  -1.50000E-01
  -5.00000E-02
 BiCGSTAB method converged.
 iter, res:           7  1.025023269341165E-011
 check the result:
i, exact solution, 1st element FEM, diff:
1   0.00000E+00   0.00000E+00   0.00000E+00
2  -2.50000E+00  -2.50000E+00   2.29092E-11
3  -4.00000E+00  -4.00000E+00   3.89715E-11
4  -4.09500E+00  -4.09500E+00   4.03242E-11
5  -4.25500E+00  -4.25500E+00   4.30216E-11
6  -4.45500E+00  -4.45500E+00   4.79270E-11
7  -4.49500E+00  -4.49500E+00   4.96323E-11
8  -4.50000E+00  -4.50000E+00   5.00506E-11

$ ./draw.py
$ display fem.png
~~~
  
TODO:
======
* do FEM by second order element(a*x^2 + b*x + c) in 1D
* try 2D, 3D

graph
======
![Alt text](fem.png?raw=true "fem")
