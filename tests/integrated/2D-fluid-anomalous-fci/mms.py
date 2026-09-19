from math import pi

from boutdata.mms import (
    DDX,
    DDZ,
    Metric,
    cos,
    diff,
    exprToStr,
    sin,
    t,
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
p = 2.0 + 0.05 * sin(4.0 * pi * x) * sin(2 * z + 0.5123)
T = p / n
D = (1.0 + 0.45 * sin(4.0 * pi * x) * cos(z + 4.231231231)) / (rho_s * rho_s * Omega_ci)
chi = (1.0 + 0.0 * x) / (rho_s * rho_s * Omega_ci)

# Turn solution into real x and z coordinates
replace = [(x, metric.x), (z, metric.z * 2.0 * pi)]

# Replace the variables with the new metric
n = n.subs(replace)
D = D.subs(replace)
p = p.subs(replace)
T = T.subs(replace)
chi = chi.subs(replace)

##############################
# Calculate time derivatives


# Density equation

dndt = (DDX(D * DDX(n)) + DDZ(D * DDZ(n))) * rho_s**2
dpdt = (
    3.0 / 2.0 * (DDX(D * T * DDX(n)) + DDZ(D * T * DDZ(n)))
    + (DDX(chi * n * DDX(T)) + DDZ(chi * n * DDZ(T)))
) * rho_s**2

#############################
# Calculate sources

Sn = diff(n, t) - dndt
Sp = diff(p, t) - 2.0 / 3.0 * dpdt
# Substitute back to get input y coordinates
replace = [(metric.x, x), (metric.z, z / (2.0 * pi))]

# Recalculate the variables
n = n.subs(replace)
p = p.subs(replace)
Sn = Sn.subs(replace)
Sp = Sp.subs(replace)
D = D.subs(replace)
chi = chi.subs(replace)
print("anomalous_D = " + exprToStr(D * (rho_s * rho_s * Omega_ci)))
print("anomalous_chi = " + exprToStr(chi * (rho_s * rho_s * Omega_ci)))

print("[Nh+]")
print("solution = " + exprToStr(n))
print("\nsource = " + exprToStr(Sn))

print("[Ph+]")
print("solution = " + exprToStr(p))
print("\nsource = " + exprToStr(Sp))
