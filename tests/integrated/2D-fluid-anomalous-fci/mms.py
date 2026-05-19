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
n = 1.0 + 0.1 * sin(2.0 * pi * x) * sin(2 * z+ 0.231)
p = 2.0 + 0.05 * sin(4.0 * pi * x) * sin(2 * z + 0.5123)

T = p / n

D = (1.0 +  0.1*sin(2.0 * pi * x) * cos(2*z))/(rho_s*rho_s*Omega_ci)
chi = (1.0 + 0.0*x) / (rho_s * rho_s * Omega_ci)




replace = [(x, metric.x), (z, metric.z * 2.0 * pi ) ]
n = n.subs(replace)
D = D.subs(replace)
chi = chi.subs(replace)
p = p.subs(replace)
T = T.subs(replace)
##############################
# Calculate time derivatives


# Density equation

dndt = (DDX(D*DDX(n)) + DDZ(D*DDZ(n))) * rho_s**2   
dpdt = ( 3.0 / 2.0 * ( DDX(D*T*DDX(n))+DDZ(D*T*DDZ(n)) ) + (DDX(chi*n*DDX(T)) + DDZ(chi*n*DDZ(T)) )  ) * rho_s**2
#############################
# Calculate sources

Sn = diff(n, t) - dndt
Sp = diff(p,t) - 2.0 / 3.0 * dpdt
# Substitute back to get input y coordinates
replace = [(metric.x, x), (metric.z, z / (2.0 * pi) ) ]
n = n.subs(replace)
p = p.subs(replace)
Sn = Sn.subs(replace)
Sp = Sp.subs(replace)
print("[n]")
print("solution = " + exprToStr(n))
print("\nsource = " + exprToStr(Sn))

print("[p]")
print("solution = " + exprToStr(p))
print("\nsource = " + exprToStr(Sp))
