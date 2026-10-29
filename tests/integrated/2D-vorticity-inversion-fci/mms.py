from math import pi

from boutdata.mms import (
    DDX,
    DDZ,
    Metric,
    exprToStr,
    sin,
    x,
    z,
)

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
Omega_ci = qe * Bnorm / (1836.0 * Me)
rho_s = 0.00022847

# Define solution in terms of input x,y,z
omega = 0.0001
n = 1.0 + 0.35 * sin(2.0 * pi * x) * sin(2 * z + 1.321312)
Te = 2.0 + 0.15 * sin(2.0 * pi * x) * sin(2 * z + 2.9911312)
Ti = 2.5 + 0.35 * sin(2.0 * pi * x) * sin(2 * z + 3.6123312)
B = 1.0 + 0.0 * x
phi = 3.35 * sin(2.0 * pi * x) * sin(4 * z + 0.8623312)

# Turn solution into real x and z coordinates
replace = [(x, metric.x), (z, metric.z * 2.0 * pi)]

# Replace the variables with the new metric
n = n.subs(replace)
Te = Te.subs(replace)
Ti = Ti.subs(replace)
B = B.subs(replace)
phi = phi.subs(replace)


def div_a_grad_perp(a, b):
    return (DDX(a * DDX(b)) + DDZ(a * DDZ(b))) * rho_s**2


vort = div_a_grad_perp(1.0 / (B * B), phi)

# Substitute back to get input y coordinates
replace = [(metric.x, x), (metric.z, z / (2.0 * pi))]

vort = vort.subs(replace)
phi = phi.subs(replace)
print("Vort = " + exprToStr(vort))
print("NVh+ = " + exprToStr(phi))
