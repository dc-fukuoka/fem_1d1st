#!/usr/bin/env python3

import f90nml
import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt

nml_file = "fort.11"
x_file   = "fort.7"
fem_file = "result.dat"

nml    = f90nml.read(nml_file)
L      = nml["params"]["length"]
lam    = nml["params"]["lambda"]
nlines = len(open(x_file).readlines())

x_coords       = np.zeros(nlines+2, dtype=np.float64)
x_coords[1:-1] = np.loadtxt(x_file, dtype=np.float64)
x_coords[-1]   = L

with open(fem_file, "rb") as fp:
    fem = np.fromfile(fp, np.float64, -1)

# exact solution: u(x) = lambda/2*x(x-2L)
exact = lam/2*x_coords*(x_coords-2*L)

sp = scipy.interpolate.InterpolatedUnivariateSpline(x_coords, exact)
sx = np.linspace(x_coords[0], x_coords[-1], 1000)
sy = sp(sx)

plt.plot(x_coords, fem, label="fem")
plt.plot(sx, sy, label="exact")
plt.xlabel("x")
plt.ylabel("y")
plt.title("1st order element FEM vs exact solution")
plt.legend(loc="upper right")
plt.grid(True)
#plt.show()
plt.savefig("fem.png")
