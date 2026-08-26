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
# n = 1 + 0.1*sin(2 * pi * x) * sin(3 * z)# * sin(t*omega)
n = 1.0 + 0.35 * sin(2.0 * pi * x) * sin(2 * z + 1.321312)
D = (1.0 + 0.45 * sin(4.0 * pi * x) * cos(z + 4.231231231)) / (rho_s * rho_s * Omega_ci)


# Turn solution into real x and z coordinates
# replace = [(x, metric.x), (z, metric.z / (2.0 * pi) ) ]
# replace = [ (y, metric.y) ]
# replace = [ (y, metric.y* Ly / (2.0 * pi)) ]


replace = [(x, metric.x), (z, metric.z * 2.0 * pi)]
n = n.subs(replace)
D = D.subs(replace)


##############################
# Calculate time derivatives


# Density equation

dndt = (DDX(D * DDX(n)) + DDZ(D * DDZ(n))) * rho_s**2

#############################
# Calculate sources

Sn = diff(n, t) - dndt

# Substitute back to get input y coordinates
replace = [(metric.x, x), (metric.z, z / (2.0 * pi))]
n = n.subs(replace)
Sn = Sn.subs(replace)
D = D.subs(replace)
print("anomalous_D = " + exprToStr(D * (rho_s * rho_s * Omega_ci)))
print("[Nh+]")
print("solution = " + exprToStr(n))
print("\nsource = " + exprToStr(Sn))
