from __future__ import division
from __future__ import print_function

from boutdata.mms import Metric, sin, cos, Div_par, Grad_par, exprToStr, diff, y, t, DDY, x, z, DDX, DDZ, Delp2, D2DX2
from math import pi

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
rho_s = 0.00022847

# Define solution in terms of input x,y,z
omega = 0.0001
#n = 1 + 0.1*sin(2 * pi * x) * sin(3 * z)# * sin(t*omega)
n = 1.0 + 0.1 * sin(2.0 * pi * x) * sin(2 * z) 
D = (1.0 + 0.0*x) / (rho_s * rho_s * Omega_ci)


# Turn solution into real x and z coordinates
#replace = [(x, metric.x), (z, metric.z / (2.0 * pi) ) ]
#replace = [ (y, metric.y) ]
#replace = [ (y, metric.y* Ly / (2.0 * pi)) ]



replace = [(x, metric.x), (z, metric.z * 2.0 * pi ) ]
n = n.subs(replace)
D = D.subs(replace)


##############################
# Calculate time derivatives


# Density equation
#dndt = (DDX(D*DDX(n)) + DDZ(D*DDZ(n))) * rho_s**2
dndt = D * Delp2(n) * rho_s**2   

#############################
# Calculate sources

Sn = diff(n, t) - dndt

# Substitute back to get input y coordinates
replace = [(metric.x, x), (metric.z, z / (2.0 * pi) ) ]
n = n.subs(replace)
Sn = Sn.subs(replace)

print("[n]")
print("solution = " + exprToStr(n))
print("\nsource = " + exprToStr(Sn))

