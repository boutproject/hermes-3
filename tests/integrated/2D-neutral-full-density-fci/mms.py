from __future__ import division
from __future__ import print_function

from boutdata.mms import Metric, sin, cos, Div_par, Grad_par, exprToStr, diff, y, t, DDY, x, z, DDX, DDZ, Delp2, D2DX2, sqrt
from math import pi
from sympy import log
# Length of the y domain
Ly = 2.0 * pi

# Atomic mass number
AA = 1.0

# metric tensor
metric = Metric()  # Identity

qe = 1.60217663e-19
Me = 9.1093837e-31
e0 = 8.85418781e-12
Pi = pi

Tnorm = 5
Nnorm = 1e18
Bnorm = 1.0
Omega_ci = qe * Bnorm / (1836.0*Me);
seconds = 1.0 / Omega_ci;
rho_s = 0.00022847
neutral_lmax = 0.02/ rho_s
# Define solution in terms of input x,y,z


Nn = (1e17 + 2e16* sin(4.0 * pi * x) * sin(2 * z + 1.33221) ) / Nnorm
Vn_x = 30 * sin(4.0 * pi * x) * sin(2 * z + 1.33221) / (rho_s/seconds)
Vn_z = 20 * sin(2.0 * pi * x) * sin(4 * z + 2.34123135) / (rho_s/seconds)

anomalous_D = 0.5 + 0.2 * sin(2.0 * pi * x) * sin(4 * z + 4.661) / ((rho_s * rho_s) / seconds)

replace = [(x, metric.x), (z, metric.z * 2.0 * pi ) ]

Nn = Nn.subs(replace)
Vn_x = Vn_x.subs(replace)
Vn_z = Vn_z.subs(replace)
anomalous_D = anomalous_D.subs(replace)



dNndt = ( -DDX(Nn*Vn_x) -DDZ(Nn*Vn_z)) * rho_s + (DDX(anomalous_D * DDX(Nn)) + DDZ(anomalous_D * DDZ(Nn))) * (rho_s * rho_s)

SNn = diff(Nn, t) - dNndt

replace = [(metric.x, x), (metric.z, z / (2.0 * pi) ) ]

Nn = Nn.subs(replace)
Vn_x = Vn_x.subs(replace)
Vn_z = Vn_z.subs(replace)
SNn = SNn.subs(replace)
anomalous_D = anomalous_D.subs(replace)

print("anomalous_D = " + exprToStr(anomalous_D * rho_s * rho_s / seconds))
print("initial_Vn_x = " + exprToStr(Vn_x * (rho_s/seconds)))
print("initial_Vn_z = " + exprToStr(Vn_z * (rho_s/seconds)))

print("[Nh]")
print("solution = " + exprToStr(Nn))
print("\nsource = " + exprToStr(SNn))
print("bndry_all = dirichlet_o2(`Nh:solution`)")
