from math import pi

from boutdata.mms import (
    DDX,
    DDZ,
    Metric,
    diff,
    exprToStr,
    sin,
    t,
    x,
    z,
)

# Length of the y domain
Ly = 2.0 * pi

# metric tensor
metric = Metric()  # Identity

qe = 1.60217663e-19
Me = 9.1093837e-31
Mp = 1836.0 * Me
e0 = 8.85418781e-12

Tnorm = 5
Nnorm = 1e18
Bnorm = 1.0
Omega_ci = qe * Bnorm / (1836.0 * Me)
rho_s = 0.00022847

lamnorm = Tnorm * Omega_ci / qe / Nnorm


def div_a_grad_perp(a, b):
    return (DDX(a * DDX(b)) + DDZ(a * DDZ(b))) * rho_s**2


# Define solution in terms of input x,y,z
phi1 = 0.35 * sin(2.0 * pi * x) * sin(2 * z + 1.321312)
Vort = 1e-2 * sin(2.0 * pi * x) * sin(2 * z - 2.521312)

B = (1.0 + 0.01 * (x - 0.5)) / Bnorm

n = 1.0 + 0.2 * sin(2.0 * pi * x) * sin(2 * z + 4.7)
T = 1.0 + 0.0 * x
AA = 1.0
Z = 1.0

viscosity = (1.0 + 0.0 * x) / (rho_s * rho_s * Omega_ci)

lam1 = 3e8 / lamnorm
lam2 = 1.0

# Turn solution into real x and z coordinates
replace = [(x, metric.x), (z, metric.z * 2.0 * pi)]

# Replace the variables with the new metric

phi1 = phi1.subs(replace)
Vort = Vort.subs(replace)

B = B.subs(replace)

n = n.subs(replace)
T = T.subs(replace)

viscosity = viscosity.subs(replace)

##############################
# Calculate some other things now
phi = phi1 / lam2
P = n * T
permittivity = e0 / ((qe * qe) * Nnorm / Mp / (Omega_ci * Omega_ci))
diamag = 0.0

##############################
# Calculate time derivatives

dVortdt = div_a_grad_perp(viscosity, Vort)

vort_from_phi = (
    div_a_grad_perp(permittivity, phi)
    + div_a_grad_perp(AA * n / (B * B), phi)
    + diamag * div_a_grad_perp(AA / (B * Z), P)
)

dphi1dt = lam1 * (vort_from_phi - Vort)

#############################
# Calculate sources

SVort = diff(Vort, t) - dVortdt

Sphi1 = diff(phi1, t) - dphi1dt


# Substitute back to get input y coordinates
replace = [(metric.x, x), (metric.z, z / (2.0 * pi))]

# Recalculate the variables
n = n.subs(replace)
T = T.subs(replace)


Vort = Vort.subs(replace)
phi1 = phi1.subs(replace)
SVort = SVort.subs(replace)
Sphi1 = Sphi1.subs(replace)


print("permittivity = " + exprToStr(permittivity))
print("density = " + exprToStr(n * Nnorm))
print("temperature = " + exprToStr(T * Tnorm))

print("[Vort]")
print("solution = " + exprToStr(Vort))
print("source = " + exprToStr(SVort))
print("bndry_all = dirichlet(Vort:solution)")
print("[phi1]")
print("solution = " + exprToStr(phi1))
print("source = " + exprToStr(Sphi1))
print("bndry_all = dirichlet(phi1:solution)")
