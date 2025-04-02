#!/usr/bin/env python3
#
# Post-process heat-conduction output, to check the energy balance

import argparse

parser = argparse.ArgumentParser(description="Analyse heat and particle fluxes")
parser.add_argument("gridfile", type=str, help="The input grid file")
parser.add_argument("datapath", type=str, help="The path to the data (BOUT.dmp files)")

args = parser.parse_args()
gridfilepath = args.gridfile
path = args.datapath

import xarray
import xhermes

yboundaries = False

bd = xhermes.open(path, geometry = "toroidal", gridfilepath = gridfilepath,
                  keep_yboundaries=yboundaries).isel(t=-1)

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



def sheath_boundary_simple(bd, gamma_e, gamma_i, Ne, Te, Ti, Zi=1, AA=1, sheath_ion_polytropic=1.0,
                           include_convective=True):
    """
    Calculate the electron and ion heat flux at the sheath, using the formula used in the
    sheath_boundary_simple component, assuming a single ion species
    with charge Zi (hydrogen=1) and atomic mass AA (hydrogen=1)

    # Returns

    flux_down, flux_up

    With units of Watts, i.e the power flowing out of each cell

    Slices at lower Y and upper Y respectively, giving heat conduction through sheath.
    Note: These do *not* include the convective heat flux, since that would usually
    be calculated in the pressure evolution (evolve_pressure component).
    """

    if not include_convective:
        gamma_e = gamma_e - 2.5
        gamma_i = gamma_i - 3.5
    
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
    q_e = ((gamma_e - 2.5) * qe * tesheath - 0.5 * me * vesheath**2) * nesheath * vesheath
    q_i = gamma_i * qe * tisheath * nesheath * vesheath

    # Multiply by cell area to get power
    flux_down_e = q_e * dx.isel(theta=y) * dz.isel(theta=y) * (J.isel(theta=y) + J.isel(theta=ym)) / (np.sqrt(g_22.isel(theta=y)) + np.sqrt(g_22.isel(theta=ym)))
    flux_down_i = q_i * dx.isel(theta=y) * dz.isel(theta=y) * (J.isel(theta=y) + J.isel(theta=ym)) / (np.sqrt(g_22.isel(theta=y)) + np.sqrt(g_22.isel(theta=ym)))

    ions_down = nesheath * vesheath * dx.isel(theta=y) * dz.isel(theta=y) * (J.isel(theta=y) + J.isel(theta=ym)) / (np.sqrt(g_22.isel(theta=y)) + np.sqrt(g_22.isel(theta=ym)))
    
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
    
    q_e = (gamma_e * qe * tesheath - 0.5 * me * vesheath**2) * nesheath * vesheath
    q_i = gamma_i * qe * tisheath * nesheath * vesheath

    flux_up_e = q_e * dx.isel(theta=y) * dz.isel(theta=y) * (J.isel(theta=y) + J.isel(theta=yp)) / (np.sqrt(g_22.isel(theta=y)) + np.sqrt(g_22.isel(theta=yp)))
    flux_up_i = q_i * dx.isel(theta=y) * dz.isel(theta=y) * (J.isel(theta=y) + J.isel(theta=yp)) / (np.sqrt(g_22.isel(theta=y)) + np.sqrt(g_22.isel(theta=yp)))

    ions_up = nesheath * vesheath * dx.isel(theta=y) * dz.isel(theta=y) * (J.isel(theta=y) + J.isel(theta=yp)) / (np.sqrt(g_22.isel(theta=y)) + np.sqrt(g_22.isel(theta=yp)))

    return flux_down_e, flux_up_e, flux_down_i, flux_up_i, ions_down, ions_up

qe = 1.602e-19 # Elementary charge [C]

# Get the conduction coefficient and fixed density from the options [m^2/s]
chi_e = bd.options["e"]["anomalous_chi"]
chi_i = bd.options["d+"]["anomalous_chi"]
D = bd.options["e"]["anomalous_D"]

Te = bd["Te"]  # Electron temperature
Ti = bd["Td+"]  # Ion temperature
Ne = bd["Ne"]  # Electron density

# Calculate radial heat fluxes at cell edges
F_L_ex, F_R_ex = Div_a_Grad_perp_upwind(bd, chi_e * Ne, qe * Te)
F_L_eD, F_R_eD = Div_a_Grad_perp_upwind(bd, qe * D * Te, Ne)
F_L_e = F_L_ex + F_L_eD
F_R_e = F_R_ex + F_R_eD

F_L_ix, F_R_ix = Div_a_Grad_perp_upwind(bd, chi_i * Ne, qe * Ti)
F_L_iD, F_R_iD = Div_a_Grad_perp_upwind(bd, qe * D * Ti, Ne)
F_L_i = F_L_ix + F_L_iD
F_R_i = F_R_ix + F_R_iD

if yboundaries:
    theta_slice = slice(2, -2)
else:
    theta_slice = slice(0, None)

# Power into domain through left boundary
# Exclude Y guard cells
total_F_L_e = F_L_e.isel(x=2, zeta=0, theta=theta_slice).sum('theta')
total_F_L_i = F_L_i.isel(x=2, zeta=0, theta=theta_slice).sum('theta')

# Power out of domain through right boundary
# Exclude Y guard cells
total_F_R_e = F_R_e.isel(x=-3, zeta=0, theta=theta_slice).sum('theta')
total_F_R_i = F_R_i.isel(x=-3, zeta=0, theta=theta_slice).sum('theta')

Pin = float(total_F_L_e + total_F_L_i)
Poutx = float(total_F_R_e + total_F_R_i)
print("Power in through X inner Pin = {} W  [{} ion {} electron]".format(
    Pin, float(total_F_L_i), float(total_F_L_e)))
print("Power out through X outer Pout,x {} W [{} ion {} electron]".format(
    Poutx, float(total_F_R_i), float(total_F_R_e)))

# Sheath

if 'sheath_boundary_simple' in bd.options:
    gamma_e = bd.options['sheath_boundary_simple']['gamma_e']
    gamma_i = bd.options['sheath_boundary_simple']['gamma_i']

    sheath_down_e, sheath_up_e, sheath_down_i, sheath_up_i, sheath_flow_down, sheath_flow_up = sheath_boundary_simple(bd, gamma_e, gamma_i, Ne, Te, Ti)
    total_sheath_down_e = sheath_down_e.isel(zeta=0, x=slice(2,-2)).sum('x')
    total_sheath_up_e = sheath_up_e.isel(zeta=0, x=slice(2,-2)).sum('x')
    total_sheath_down_i = sheath_down_i.isel(zeta=0, x=slice(2,-2)).sum('x')
    total_sheath_up_i = sheath_up_i.isel(zeta=0, x=slice(2,-2)).sum('x')

    Poutd = float(total_sheath_down_e + total_sheath_down_i)
    Poutu = float(total_sheath_up_e + total_sheath_up_i)
    print("Power out through lower sheath: Pout,d = {} W [{} ion {} electron]".format(
        Poutd, float(total_sheath_down_i), float(total_sheath_down_e)))
    print("Power out through upper sheath: Pout,u = {} W [{} ion {} electron]".format(
        Poutu, float(total_sheath_up_i), float(total_sheath_up_e)))
else:
    Poutd = 0.0
    Poutu = 0.0

# Energy content

cell_volume = bd['J'] * bd['dx'] * bd['dy'] * bd['dz']
domain_volume = float(cell_volume.isel(x=slice(2,-2), theta=theta_slice).sum())  # In m^3, excluding X guard cells

print("Domain volume: {} m^3".format(domain_volume))

# Atomic processes, radiation and other sources/sinks of thermal energy

Rex = bd['Rd+_ex']
Rrec = bd['Rd+_rec']

excitation_radiation = 1.5 * (Rex * cell_volume).isel(x=slice(2,-2), theta=theta_slice).sum(['x', 'theta', 'zeta'])  # In Watts
recombination_radiation = 1.5 * (Rrec * cell_volume).isel(x=slice(2,-2), theta=theta_slice).sum(['x', 'theta', 'zeta'])  # In Watts

Prad = -float(excitation_radiation + recombination_radiation)
print("Net power radiation: Prad = {} W [{} excitation {} recombination]".format(
    Prad, -float(excitation_radiation), -float(recombination_radiation)))

Pnet = Pin - Poutx - Poutd - Poutu - Prad
print("Net input power Pnet = Pin - Pout,x - Pout,d - Pout,u - Prad = {} W".format(Pnet))

print("------------")

F_L, F_R = Div_a_Grad_perp_upwind(bd, D * Ne / Ne, Ne)
total_particles_in_L = float(F_L.isel(x=2, zeta=0, theta=theta_slice).sum('theta'))
total_particles_out_R = float(F_R.isel(x=-3, zeta=0, theta=theta_slice).sum('theta'))

sheath_flow_down = float(sheath_flow_down.isel(zeta=0, x=slice(2,-2)).sum('x'))
sheath_flow_up = float(sheath_flow_up.isel(zeta=0, x=slice(2,-2)).sum('x'))

print("Particle flux in through X inner: {} /s".format(total_particles_in_L))
print("Particle flux out through X outer: {} /s".format(total_particles_out_R))

print("Particle flux to lower sheath: {} /s".format(sheath_flow_down))
print("Particle flux to upper sheath: {} /s".format(sheath_flow_up))

frecyc = 1. - (total_particles_in_L - total_particles_out_R) / (sheath_flow_down + sheath_flow_up)
print("Recycling fraction: {}".format(frecyc))

Siz = bd['Sd+_iz']
Srec = bd['Sd+_rec']

total_iz = float((Siz * cell_volume).isel(x=slice(2,-2), theta=theta_slice).sum(['x', 'theta', 'zeta']))
total_rec = float((Srec * cell_volume).isel(x=slice(2,-2), theta=theta_slice).sum(['x', 'theta', 'zeta']))

print("Total ionization rate: {} /s".format(total_iz))
print("Total recombination rate: {} /s".format(total_rec))

print("------------")

# Calculate rate of change of energy
ddt_Ep = 0.0
ddt_Ek = 0.0
for species in ['e', 'd+', 'd']:
    # Internal energy
    ddt_Ep += 1.5 * (bd['ddt(P' + species + ')'] * cell_volume).isel(x=slice(2,-2), theta=theta_slice).sum(['x', 'theta', 'zeta'])
    # Kinetic energy
    ddt_Ek += bd['V' + species] * bd['ddt(NV' + species + ')'] - 0.5 * bd['NV' + species] * bd['V'+species] * bd['ddt(N'+species+')'] / bd['N'+species]
ddt_Ep = float(ddt_Ep)
ddt_Ek = float(ddt_Ek)

print("Rate of change of energy: {} [{} internal {} kinetic]".format(ddt_Ep + ddt_Ek, ddt_Ep, ddt_Ek))
