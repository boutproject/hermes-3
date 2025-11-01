from __future__ import division
from __future__ import print_function
from sympy import log
from boutdata.mms import Metric, sin, cos, Div_par, Grad_par, exprToStr, diff, y, t, DDY, sqrt
from math import pi
import numpy as np
# Length of the y domain
Ly = 2.0 * pi

# Atomic mass number
AA = 1.0 / 1836.0

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
print(f"Omega_ci = {np.round(Omega_ci,4)}")
print(f"seconds = {1.0/Omega_ci}")
# Define solution in terms of input x,y,z
omega = 0.0001
n = 1 + 0.3*sin(y)# * sin(t*omega)
T = 4 + 2.0*sin(y)
p = n * T

# Turn solution into real x and z coordinates
replace = [ (y, metric.y*2*pi/Ly) ]
#replace = [ (y, metric.y) ]
#replace = [ (y, metric.y* Ly / (2.0 * pi)) ]

rho_s = 0.0002284697436697996

n = n.subs(replace)
T = T.subs(replace)
p = p.subs(replace)

##############################
# Calculate time derivatives


# Density equation
dndt = 0.0
# Pressure equation
def coulog(logden, logtemp):
    return 30.4 - 0.5 * logden + (5. / 4) * logtemp - sqrt(1e-5 + (logtemp - 2.0)**2 / 16.)

def calc_frequency(den,temp):
    logden = log(den)
    logtemp = log(temp)

    
    coulomb_log = coulog(logden, logtemp)
    v1sq = 2.0 * temp * qe / Me

    nu = (qe**4) * den * coulomb_log * 2.0 / (3.0 * (Pi * 2.0 * v1sq)**1.5  * (e0*Me)**2)
    
    return  nu/Omega_ci

tau = 1.0 / calc_frequency(n * Nnorm, T * Tnorm)
kappa = (3.16/sqrt(2)) * p * tau / AA
dpdt = (2.0/3.0) * Div_par(kappa * Grad_par(T)) * (rho_s**2)

#############################
# Calculate sources

Sn = diff(n, t)# - dndt
Sp = diff(p, t) - dpdt

# Substitute back to get input y coordinates
replace = [ (metric.y, y*Ly/(2*pi) ) ]
#replace = [ (metric.y, y ) ]
#replace = [ (metric.y, y*(2*pi) / Ly ) ]

n = n.subs(replace)
p = p.subs(replace)

Sn = Sn.subs(replace)
Sp = Sp.subs(replace)


print("[n]")
print("solution = " + exprToStr(n))
print("\nsource = " + exprToStr(Sn))

print("\n[p]")
print("solution = " + exprToStr(p))
print("\nsource = " + exprToStr(Sp))
