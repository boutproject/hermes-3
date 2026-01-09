#!/usr/bin/env python3
#
# Check the poloidal integral of the radial magnetic drift, Curl(b/B)
#
# The sum over the core poloidal points of J * Curlb_Bx should be zero

import sys

if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} gridfile.nc")
    print("Note: The input grid file is modified")
    sys.exit(1)

filename = sys.argv[1]

from boututils.datafile import DataFile
import numpy as np

# Read values from current grid file
with DataFile(filename) as g:
    J = g["J"]
    B = g["Bxy"]

    # Curvature components
    bxcvx = g["bxcvx"]
    bxcvz = g["bxcvz"]

    try:
      sinty = g["sinty"]  # Integrated magnetic shear
    except KeyError:
      sinty = 0.0

    # Curl(b / B)
    curlb_Bx = 2 * bxcvx / B

    # Poloidal branch-cut indices
    j11 = g["jyseps1_1"]
    j12 = g["jyseps1_2"]
    j21 = g["jyseps2_1"]
    j22 = g["jyseps2_2"]

    # Get the separatrix index
    ix1 = g["ixseps1"]
    ix2 = g["ixseps2"]
    ixcore = np.min([ix1, ix2])

    # The coordinate system depends on the transform
    try:
        parallel_transform = g["parallel_transform"]
    except KeyError:
        parallel_transform = "identity"

# Calculate J * Curl(b/B) in the core region
Jc = J * curlb_Bx
Jc_core = np.concatenate(
    (Jc[:, (j11 + 1) : (j21 + 1)], Jc[:, (j12 + 1) : (j22 + 1)]), axis=1
)

# Sum over y (poloidal) index in the core
Jc_sum = np.sum(Jc_core, axis=1)

# Rescale curvature
# Keep places where curvature = 0, so multiply positive and negative values
# by factors.
# Positive values are multiplied by "factor"; negative values are divided
#
# Curvature = positive_curvature * factor + negative_curvature / factor
#


# Sum over just the positive part of the curvature
Jc_positive_sum = np.sum(np.clip(Jc_core, 0.0, None), axis=1)

factor = np.sqrt(1 - Jc_sum / Jc_positive_sum)

# Keep the factor constant outside the separatrix
# so that curvature is still continuous
factor_edge = factor[ixcore - 1]  # Last closed flux surface

bxcvx2 = (
    np.clip(bxcvx, 0.0, None) * factor_edge
    + np.clip(bxcvx, None, 0.0) / factor_edge  # Positive parts
)  # Negative parts

# Go through core flux surfaces
for i in range(ixcore):
    # Inner core
    bxcvx2[i, (j11 + 1) : (j21 + 1)] = (
        np.clip(bxcvx[i, (j11 + 1) : (j21 + 1)], 0.0, None) * factor[i]
        + np.clip(bxcvx[i, (j11 + 1) : (j21 + 1)], None, 0.0) / factor[i]
    )

    # Outer core
    bxcvx2[i, (j12 + 1) : (j22 + 1)] = (
        np.clip(bxcvx[i, (j12 + 1) : (j22 + 1)], 0.0, None) * factor[i]
        + np.clip(bxcvx[i, (j12 + 1) : (j22 + 1)], None, 0.0) / factor[i]
    )


# Check that the sum of the curvature sums to zero over poloidal angle

Jc2 = J * 2 * bxcvx2 / B
Jc2_core = np.concatenate(
    (Jc2[:, (j11 + 1) : (j21 + 1)], Jc2[:, (j12 + 1) : (j22 + 1)]), axis=1
)
Jc2_sum = np.sum(Jc2_core, axis=1)

import matplotlib.pyplot as plt

plt.plot(Jc_sum, label="Input ")
plt.plot(Jc2_sum, label="Adjusted")
plt.xlabel("Radial index")
plt.ylabel(r"$\sum_{y,core} J \nabla\times(\mathbf{b}/B)$")
plt.legend()
plt.show()

with DataFile(filename, write=True) as outf:
    # Update curvature
    outf["bxcvx"] = bxcvx2

    if parallel_transform == "identity":
        # Need to correct the z component
        outf["bxcvz"] = bxcvz + sinty * (bxcvx - bxcvx2)
