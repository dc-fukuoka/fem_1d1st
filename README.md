one dimensional first element finite element method(Galerkin method) using BiCGSTAB method.

~~~
$ make
$ ./a.out
 e%nelements:           3
 e%nnodes:           4
 points:
   0.00000E+00   1.00000E+00   2.00000E+00   3.00000E+00
 lhs A:
   1.00000E+00   0.00000E+00   0.00000E+00   0.00000E+00
  -1.00000E+00   2.00000E+00  -1.00000E+00   0.00000E+00
   0.00000E+00  -1.00000E+00   2.00000E+00  -1.00000E+00
   0.00000E+00   0.00000E+00  -1.00000E+00   1.00000E+00

 rhs b:
   0.00000E+00
  -1.00000E+00
  -1.00000E+00
  -5.00000E-01
 BiCGSTAB method converged.
 iter, res:           3  8.881784197001252E-016
i, exact solution, 1st element FEM, diff:           1   0.00000E+00   0.00000E+00   0.00000E+00
i, exact solution, 1st element FEM, diff:           2  -2.50000E+00  -2.50000E+00  -4.44089E-16
i, exact solution, 1st element FEM, diff:           3  -4.00000E+00  -4.00000E+00  -8.88178E-16
i, exact solution, 1st element FEM, diff:           4  -4.50000E+00  -4.50000E+00  -8.88178E-16
~~~
  
TODO:
======
* add a graph to compare the first element(patching linear functions) and exact one
* do FEM by second order element(a*x^2 + b*x + c) in 1D
* try 2D, 3D
