from __future__ import division
from __future__ import print_function

from boutdata.mms import Metric, sin, cos, Div_par, Grad_par, exprToStr, diff, y, t, DDY, sqrt
from math import pi
from sympy import log
import numpy as np
# Length of the y domain
Ly = 2.0 * pi

# Atomic mass number
AA = 1.0
charge = 1.0
qe = 1.60217663e-19
Me = 9.1093837e-31
Mp = Me * 1836.0
e0 = 8.85418781e-12
Pi = pi

Tnorm = 5
Nnorm = 1e18
Bnorm = 1.0
Omega_ci = qe * Bnorm / (1836.0*Me);
print(f"Omega_ci = {np.round(Omega_ci,4)}")
print(f"seconds = {1.0/Omega_ci}")

# metric tensor
metric = Metric()  # Identity

# Define solution in terms of input x,y,z

omega = 1e-4

n = 1 + 0.1*sin(2*y)# * sin(t*omega)
p = 1 + 0.1*cos(3*y)# * sin(t*omega)
mnv = AA * 0.1*sin(y)# * sin(2*t*omega)
B = 1.0 + 0.0 * y

# Turn solution into real x and z coordinates
replace = [ (y, metric.y*2*pi/Ly) ]
#replace = [ (y, metric.y) ]
#replace = [ (y, metric.y* Ly / (2.0 * pi)) ]

rho_s = 0.00022847

n = n.subs(replace)
p = p.subs(replace)
mnv = mnv.subs(replace)
B = B.subs(replace)

##############################
# Calculate time derivatives

nv = mnv / AA
v = nv / n
gamma = 5./3
T = p/n


def calc_frequency(den,temp):
    logden = log(den)
    logtemp = log(temp)

    charge1 = charge*qe
    density1 = den
    temp1 = temp
    mass1 = AA * Mp

    #29.91 - log((Z1 * Z2 * (AA1 + AA2)) / (AA1 * Tlim2 + AA2 * Tlim1) * sqrt(Nlim1 * SQ(Z1) / Tlim1 + Nlim2 * SQ(Z2) / Tlim2));

    coulomb_log = 29.91 - log( ( charge * charge * (AA+AA)) / (AA*temp1 + AA*temp1) * sqrt(2.0* den * charge**2 / temp))                               
    v1sq = 2.0 * temp * qe / mass1

    # SQ(charge1 * charge2) * Nlim2 * floor(coulomb_log, 1.0) * (1. + mass1 / mass2) / (3 * pow(PI * (v1sq + v2sq), 1.5) * SQ(SI::e0 * mass1));
    nu = ((charge1*charge1)**2 ) * den * coulomb_log * (1.0 + mass1/mass1) / (3.0 * (pi * 2.0 * v1sq)**1.5 * (e0*mass1)**2)                          
    return  nu/Omega_ci

tau = 1.0 / calc_frequency(n * Nnorm, T * Tnorm)

eta = 1.28 * p * tau

div_Pi_cipar = sqrt(B) * Div_par( (eta / B) * Grad_par(sqrt(B) * v) ) *	rho_s**2


# Density equation
dndt = - Div_par(nv) * rho_s

# Pressure equation
dpdt = - Div_par(p*v) * rho_s - (gamma-1.0)*p*Div_par(v) * rho_s - 2.0/3.0 * v * div_Pi_cipar

# Momentum equation
dmnvdt = - Div_par(mnv*v) * rho_s - Grad_par(p) * rho_s + div_Pi_cipar


#############################
# Calculate sources

Sn = diff(n, t) - dndt
Sp = diff(p, t) - dpdt
Smnv = diff(mnv, t) - dmnvdt

# Substitute back to get input y coordinates
replace = [ (metric.y, y*Ly/(2*pi) ) ]
#replace = [ (metric.y, y ) ]
#replace = [ (metric.y, y*(2*pi) / Ly ) ]

n = n.subs(replace)
p = p.subs(replace)
mnv = mnv.subs(replace)

Sn = Sn.subs(replace)
Sp = Sp.subs(replace)
Smnv = Smnv.subs(replace)

print("[n]")
print("solution = " + exprToStr(n))
print("\nsource = " + exprToStr(Sn))

print("\n[p]")
print("solution = " + exprToStr(p))
print("\nsource = " + exprToStr(Sp))

print("\n[nv]")
print("solution = " + exprToStr(mnv))
print("\nsource = " + exprToStr(Smnv))
