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
rho_s = 0.00022847
neutral_lmax = 0.02/ rho_s
# Define solution in terms of input x,y,z


Nn = (1e17 + 2e16* sin(4.0 * pi * x) * sin(2 * z + 1.33221) ) / Nnorm
Tn = (10.0 + 2.0 * sin(2.0 * pi * x) * sin(2 * z - 0.315312)) / Tnorm
#Tn = (10.0 + 0.0 * x ) /Tnorm
Pn = Nn * Tn
Rnn = sqrt(Tn / AA) / neutral_lmax
Dnn = (Tn / AA) / Rnn
DnnNn = Dnn * Nn
DnnPn = Dnn * Pn
kappa_n = (5.0 / 2.0) * DnnNn 
logPn = log(Pn)

replace = [(x, metric.x), (z, metric.z * 2.0 * pi ) ]

Nn = Nn.subs(replace)
Tn = Tn.subs(replace)
Pn = Pn.subs(replace)
Rnn = Rnn.subs(replace)
Dnn = Dnn.subs(replace)
DnnNn = DnnNn.subs(replace)
DnnPn = DnnPn.subs(replace)
kappa_n = kappa_n.subs(replace)
logPn = logPn.subs(replace)
##############################
# Calculate time derivatives


# Density equation

dNndt = ( DDX(DnnNn * DDX(logPn)) + DDZ(DnnNn * DDZ(logPn)) ) * rho_s**2

dPndt = (5.0 / 3.0) * ( DDX(DnnPn * DDX(logPn)) + DDZ(DnnPn * DDZ(logPn)) ) * rho_s**2

SNn = diff(Nn, t) - dNndt
SPn = diff(Pn, t) - dPndt
# Substitute back to get input y coordinates
replace = [(metric.x, x), (metric.z, z / (2.0 * pi) ) ]

Nn = Nn.subs(replace)
Pn = Pn.subs(replace)
SNn = SNn.subs(replace)
SPn = SPn.subs(replace)

print("[Nn]")
print("solution = " + exprToStr(Nn))
print("\nsource = " + exprToStr(SNn))

print("[Pn]")
print("solution = " + exprToStr(Pn))
print("\nsource = " + exprToStr(SPn))

