from __future__ import division
from __future__ import print_function

from boutdata.mms import Metric, sin, cos, Div_par, Grad_par, exprToStr, diff, y, t, DDY
from math import pi

# Length of the y domain
Ly = 2.0 * pi

# Atomic mass number
AA = 1.0

# metric tensor
metric = Metric()  # Identity

# Define solution in terms of input x,y,z
omega = 0.0001
n = 1 + 0.1*sin(2*y)# * sin(t*omega)
p = 1 + 0.1*cos(3*y)# * sin(t*omega)
mnv = AA * 0.1*sin(y)# * sin(2*t*omega)

# Turn solution into real x and z coordinates
replace = [ (y, metric.y*2*pi/Ly) ]
#replace = [ (y, metric.y) ]
#replace = [ (y, metric.y* Ly / (2.0 * pi)) ]

rho_s = 0.00022847

n = n.subs(replace)
p = p.subs(replace)
mnv = mnv.subs(replace)

##############################
# Calculate time derivatives

nv = mnv / AA
v = nv / n
gamma = 5./3

# Density equation
dndt = - Div_par(nv) * rho_s

# Pressure equation
dpdt = - Div_par(p*v) * rho_s - (gamma-1.0)*p*Div_par(v) * rho_s

# Momentum equation
dmnvdt = - Div_par(mnv*v) * rho_s - Grad_par(p) * rho_s


#############################
# Calculate sources

Sn = diff(n, t)# - dndt
Sp = diff(p, t)# - dpdt
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
