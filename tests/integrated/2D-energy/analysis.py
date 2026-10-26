import xhermes

path = "data-Eion"

bd = xhermes.open(path)

cell_volume = bd["J"] * bd["dx"] * bd["dy"] * bd["dz"]
spatial_dims = [d for d in bd.dims if d != "t"]


def volume_integrate(value):
    return (value * cell_volume).sum(dim=spatial_dims)


if "Eion" in bd:
    total_Eion = volume_integrate(bd["Eion"])

Eelectron = (3.0 / 2) * bd["Pe"]
total_Eelectron = volume_integrate(Eelectron)

for term in [
    "DivEExB",
    "DivEReynolds",
    "DivEpol",
    "phi_Div_Je",
    "Eion_src",
    "Eion",
    "ion_KE_perp",
    "ddt(Eion)",
]:
    if term in bd:
        total = volume_integrate(bd[term]).values
        print(f"{term} =  {total}")

for term in ["Pe", "ddt(Pe)"]:
    if term in bd:
        total = (3.0 / 2) * volume_integrate(bd[term]).values
        print(f"3/2 * {term} =  {total}")


ddt_energy = volume_integrate(bd["ddt(Eion)"]) + (3.0 / 2) * volume_integrate(
    bd["ddt(Pe)"]
)
print(f"ddt(Eion + Eelectron) = {ddt_energy.values}")
