import tcvx21
from tcvx21 import Quantity
import matplotlib.pyplot as plt
from tcvx21.record_c import Record
import numpy as np
from scipy.stats import linregress
from scipy.optimize import curve_fit

plt.style.use(tcvx21.style_sheet)

tcv = Record("TCV-X21/1.experimental_data/TCV_forward_field.nc", color="C0")

hermes = Record("hermes_tcvx21_data.nc", color="C1")
hermes_highres = Record("hermes_tcvx21_highres.nc", color="C2")

# Get Te, ne profiles

ne_tcv_ts = tcv.get_observable("TS", "density")
ne_tcv_fhrp = tcv.get_observable("FHRP", "density")
ne_lowres = hermes.get_observable("FHRP", "density")
ne_highres = hermes_highres.get_observable("FHRP", "density")


def make_mask(positions):
    return np.logical_and(
        positions > Quantity(0.0, "cm"),
        positions < Quantity(1.5, "cm"),
    )


def length_scale(dataset):
    mask = make_mask(dataset.positions)
    x = dataset.positions[mask].magnitude
    y = dataset.values[mask].magnitude
    y_err = dataset.errors[mask].magnitude

    # Get initial guess
    result = linregress(x, np.log(y))
    L = -1.0 / result.slope
    y0 = np.exp(result.intercept)

    # Use curve_fit to calculate errors
    def logfit(x, y0, L):
        return y0 * np.exp(-x / L)

    result = curve_fit(logfit, x, y, p0=(y0, L), sigma=y_err, absolute_sigma=True)

    y0 = result[0][0]
    L = result[0][1]
    y0_err, L_err = np.sqrt(np.diag(result[1]))

    return (
        (
            Quantity(L, dataset.positions.units),
            Quantity(L_err, dataset.positions.units),
        ),
        (Quantity(y0, dataset.values.units), Quantity(y0_err, dataset.values.units)),
    )


for label, dataset in [
    ("lambda_n (TS)", ne_tcv_ts),
    ("lambda_n (FHRP)", ne_tcv_fhrp),
    ("lambda_n (Hermes-3 lowres)", ne_lowres),
    ("lambda_n (Hermes-3 highres)", ne_highres),
    ("lambda_T (TS)", tcv.get_observable("TS", "electron_temp")),
    ("lambda_T (FHRP)", tcv.get_observable("FHRP", "electron_temp")),
    ("lambda_T (Hermes-3 lowres)", hermes.get_observable("FHRP", "electron_temp")),
    (
        "lambda_T (Hermes-3 highres)",
        hermes_highres.get_observable("FHRP", "electron_temp"),
    ),
]:
    result = length_scale(dataset)
    L, L_err = result[0]
    y0, y0_err = result[1]
    print(f"{label}: {L} +/- {L_err}\n\tIntercept: {y0} +/- {y0_err}")

# Make plots of midplane and target profiles
for region, measurement in [
    ("LFS-LP", "density"),
    ("LFS-LP", "electron_temp"),
    ("LFS-LP", "potential"),
    ("LFS-LP", "current"),
    ("LFS-LP", "vfloat"),
    ("FHRP", "density"),
    ("FHRP", "electron_temp"),
    ("FHRP", "potential"),
    ("HFS-LP", "density"),
    ("HFS-LP", "electron_temp"),
    ("HFS-LP", "potential"),
    ("HFS-LP", "current"),
    ("HFS-LP", "vfloat"),
]:
    fig, ax = plt.subplots()
    tcv.get_observable(region, measurement).plot(ax=ax)
    hermes.get_observable(region, measurement).plot(ax=ax)
    hermes_highres.get_observable(region, measurement).plot(ax=ax)
    ax.set_title(f"{region}:{measurement}")
    fig.legend()
    plt.show()
