# Post-process heat-conduction output, to check the energy balance

import xarray
import xhermes

yboundaries = False

bd = xhermes.open(".", geometry = "toroidal", gridfilepath = "../tokamak.nc",
                  keep_yboundaries=yboundaries)

import numpy as np

# Radial fluxes due to cross-field diffusion

def Div_a_Grad_perp_upwind(bd, a, f):
    """
    # Returns

    Tuple of two quantities:
    (F_L, F_R)

    These are the flow into the cell from the left (x-1),
    and the flow out of the cell to the right (x+1). 

    Note: *NOT* the flux; these are already integrated over cell boundaries

    # Example

    F_L, F_R = Div_a_Grad_perp_upwind(chi * N, e * T)

    - chi is the heat diffusion in m^2/2
    - N is density in m^-3
    - T is temperature in eV

    Then F_L and F_R would have units of Watts,

    The time-derivative of the pressure in that cell would be:

    d/dt (3/2 P) = (F_L - F_R) / V

    where V = dx * dy * dz * J is the volume of the cell
    
    """

    J = bd["J"] # Jacobian
    g11 = bd["g11"]
    dx = bd["dx"]
    dy = bd["dy"]
    dz = bd["dz"]

    F_R = xarray.zeros_like(f) # Flux to the right
    F_L = xarray.zeros_like(f) # Flux from the left

    for x in bd.x[:-1]:
        xp = x + 1  # The next X cell
        # Note: Order of array operations matters for shape of the result
        gradient = (f.isel(x=xp) - f.isel(x=x)) * (J.isel(x=x) * g11.isel(x=x) + J.isel(x=xp) * g11.isel(x=xp))  / (dx.isel(x=x) + dx.isel(x=xp))

        flux = -gradient * 0.5*(a.isel(x=x) + a.isel(x=xp))
        
        # if gradient > 0:
        #     # Flow from x+1 to x
        #     flux = -gradient * a.isel(x=xp)  # Note: Negative flux = flow into this cell from right
        # else:
        #     # Flow from x to x+1
        #     flux = -gradient * a.isel(x=x)  # Positive flux => Flow from this cell to the right

        # Need to multiply by dy * dz because these are assumed constant in X in the calculation
        # of flux and cell volume.
        flux *= dy.isel(x=x) * dz.isel(x=x)

        F_R[dict(x=x)] = flux
        F_L[dict(x=xp)] = flux

    return F_L, F_R



def sheath_boundary_simple(bd, gamma_e, Ne, Te, Ti, Zi=1, AA=1, sheath_ion_polytropic=1.0):
    """
    Calculate the electron heat flux at the sheath, using the formula used in the
    sheath_boundary_simple component, assuming a single ion species
    with charge Zi (hydrogen=1) and atomic mass AA (hydrogen=1)

    # Returns

    flux_down, flux_up

    With units of Watts, i.e the power flowing out of each cell

    Slices at lower Y and upper Y respectively, giving heat conduction through sheath.
    Note: These do *not* include the convective heat flux, since that would usually
    be calculated in the pressure evolution (evolve_pressure component).
    """

    J = bd['J']
    dx = bd['dx']
    dy = bd['dy']
    dz = bd['dz']
    g_22 = bd['g_22']

    # Lower y
    if yboundaries:
        y = 2 # First point in the domain
        ym = y - 1
        Ne_m = Ne.isel(theta=ym)
        Te_m = Te.isel(theta=ym)
        Ti_m = Ti.isel(theta=ym)
    else:
        y = 0
        ym = y # Same metric tensor component in boundary cells as in domain
        yp = y + 1 # For extrapolating boundary
        
        Ne_m = Ne.isel(theta=y)**2 / Ne.isel(theta=yp)
        Te_m = Te.isel(theta=y)**2 / Te.isel(theta=yp)
        Ti_m = Ti.isel(theta=y)**2 / Ti.isel(theta=yp)

    nesheath = 0.5 * (Ne.isel(theta=y) + Ne_m)
    tesheath = 0.5 * (Te.isel(theta=y) + Te_m)
    tisheath = 0.5 * (Ti.isel(theta=y) + Ti_m)

    qe = 1.602e-19 # Elementary charge [C]
    mp = 1.67e-27 # Proton mass [kg]
    me = 9.11e-31 # Electron mass [kg]
    
    # Ion flow speed
    C_i = np.sqrt((sheath_ion_polytropic * qe * tisheath + Zi * qe * tesheath) / (AA * mp))

    vesheath = C_i  # Assuming no current

    # Parallel heat flux in W/m^2.
    # Note: Corrected for 5/2Pe convective thermal flux, and small electron kinetic energy flux
    # so gamma_e is the total *energy* flux coefficient.
    q = ((gamma_e - 2.5) * qe * tesheath - 0.5 * me * vesheath**2) * nesheath * vesheath

    # Multiply by cell area to get power
    flux_down = q * dx.isel(theta=y) * dz.isel(theta=y) * (J.isel(theta=y) + J.isel(theta=ym)) / (np.sqrt(g_22.isel(theta=y)) + np.sqrt(g_22.isel(theta=ym)))

    # Repeat for upper Y boundary
    if yboundaries:
        y = -3 # First point in the domain
        yp = y + 1
        Ne_p = Ne.isel(theta=yp)
        Te_p = Te.isel(theta=yp)
        Ti_p = Ti.isel(theta=yp)
    else:
        y = -1
        yp = y # Same metric tensor component in boundary cells as in domain
        ym = y - 1 # For extrapolating boundary
        
        Ne_p = Ne.isel(theta=y)**2 / Ne.isel(theta=ym)
        Te_p = Te.isel(theta=y)**2 / Te.isel(theta=ym)
        Ti_p = Ti.isel(theta=y)**2 / Ti.isel(theta=ym)

    nesheath = 0.5 * (Ne.isel(theta=y) + Ne_p)
    tesheath = 0.5 * (Te.isel(theta=y) + Te_p)
    tisheath = 0.5 * (Ti.isel(theta=y) + Ti_p)
    C_i = np.sqrt((sheath_ion_polytropic * qe * tisheath + Zi * qe * tesheath) / (AA * mp))
    vesheath = C_i
    q = ((gamma_e - 2.5) * qe * tesheath - 0.5 * me * vesheath**2) * nesheath * vesheath
    flux_up = q * dx.isel(theta=y) * dz.isel(theta=y) * (J.isel(theta=y) + J.isel(theta=yp)) / (np.sqrt(g_22.isel(theta=y)) + np.sqrt(g_22.isel(theta=yp)))

    return flux_down, flux_up

qe = 1.602e-19 # Elementary charge [C]

# Get the conduction coefficient and fixed density from the options [m^2/s]
chi = bd.options["e"]["anomalous_chi"]

Te = bd["Te"]  # Electron temperature
Ne = xarray.full_like(Te, bd.options["e"]["density"]) # m^-3

# Calculate radial heat fluxes at cell edges
F_L, F_R = Div_a_Grad_perp_upwind(bd, chi * Ne, qe * Te)

if yboundaries:
    theta_slice = slice(2, -2)
else:
    theta_slice = slice(0, None)

# Power into domain through left boundary
# Exclude Y guard cells
total_F_L = F_L.isel(x=2, zeta=0, theta=theta_slice).sum('theta')

# Power out of domain through right boundary
# Exclude Y guard cells
total_F_R = F_R.isel(x=-3, zeta=0, theta=theta_slice).sum('theta')

Pin = float(total_F_L[-1])
Poutx = float(total_F_R[-1])
print("Power in through X inner Pin = {} W".format(Pin))
print("Power out through X outer Pout,x {} W".format(Poutx))

# Sheath

if 'sheath_boundary_simple' in bd.options:
    gamma_e = bd.options['sheath_boundary_simple']['gamma_e']
    Ti = Te

    sheath_down, sheath_up = sheath_boundary_simple(bd, gamma_e, Ne, Te, Ti)
    total_sheath_down = sheath_down.isel(zeta=0, x=slice(2,-2)).sum('x')
    total_sheath_up = sheath_up.isel(zeta=0, x=slice(2,-2)).sum('x')

    Poutd = float(total_sheath_down.isel(t=-1))
    Poutu = float(total_sheath_up.isel(t=-1))
    print("Power out through lower sheath: Pout,d = {} W".format(Poutd))
    print("Power out through upper sheath: Pout,u = {} W".format(Poutu))
else:
    Poutd = 0.0
    Poutu = 0.0

print("------------")
Pnet = Pin - Poutx - Poutd - Poutu
print("Net input power Pnet = Pin - Pout,x - Pout,d - Pout,u = {} W".format(Pnet))

# Energy content

cell_volume = bd['J'] * bd['dx'] * bd['dy'] * bd['dz']
domain_volume = float(cell_volume.isel(x=slice(2,-2), theta=theta_slice).sum())  # In m^3, excluding X guard cells

Pe = bd["Pe"]
energy_content = 1.5 * (Pe * cell_volume).isel(x=slice(2,-2), theta=theta_slice).sum(['x', 'theta', 'zeta'])  # In Joules
dWdt = energy_content.differentiate('t')

print("Volume of domain: {} m^3".format(domain_volume))
print("Rate of change of energy content dW/dt = {} W".format(float(dWdt[-1])))

print("------------")
print("Power imbalance Pnet - dW/dt = {} W".format(Pnet - float(dWdt[-1])))
