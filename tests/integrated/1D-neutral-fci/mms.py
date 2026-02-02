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


Nn = (1e17 + 2e16* sin(y + 1.74773) ) / Nnorm
Tn = (10.0 + 2.0 * sin(2.0 * y + 5.23131)) / Tnorm
Vn = (0.0 + 0.1 * sin(y))
NVn = AA * Nn * Vn
#Tn = (10.0 + 0.0 * x ) /Tnorm
Pn = Nn * Tn
Rnn = sqrt(Tn / AA) / neutral_lmax
Dnn = (Tn / AA) / Rnn
DnnNn = Dnn * Nn
DnnPn = Dnn * Pn
kappa_n = (5.0 / 2.0) * DnnNn 
logPn = log(Pn)

replace = [ (y, metric.y*2*pi/Ly) ]

Nn = Nn.subs(replace)
Tn = Tn.subs(replace)
Vn = Vn.subs(replace)
NVn = NVn.subs(replace)
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



dNndt = -Div_par(Nn * Vn) * rho_s
DNVndt = (-AA * Div_par(NVn * Vn) - Grad_par(Pn)) * rho_s
dPndt = (-(5.0 / 3.0) * Div_par(Pn * Vn) - (2.0 / 3.0) * Vn * Grad_par(Pn)) * rho_s



SNn = diff(Nn, t) - dNndt
SNVn = diff(NVn, t) - DNVndt
SPn = diff(Pn, t) - dPndt



replace = [ (metric.y, y*Ly/(2*pi) ) ]



Nn = Nn.subs(replace)
NVn = NVn.subs(replace)
Pn = Pn.subs(replace)
SNn = SNn.subs(replace)
SNVn = SNVn.subs(replace)
SPn = SPn.subs(replace)

print("[Nn]")
print("solution = " + exprToStr(Nn))
print("\nsource = " + exprToStr(SNn))

print("[NVn]")
print("solution = " + exprToStr(NVn))
print("\nsource = " + exprToStr(SNVn))

print("[Pn]")
print("solution = " + exprToStr(Pn))
print("\nsource = " + exprToStr(SPn))

