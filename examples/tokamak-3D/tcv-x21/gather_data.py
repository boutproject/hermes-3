import numpy as np
from boutdata import collect
from boututils.datafile import DataFile
import pickle

qe = 1.602e-19
AA_me = 1.0 / 1836  # Electron mass divided by proton mass


def extract_data(path, gridfilepath, ymid=18):
    """
    Read data from a Hermes-3/BOUT++ output directory.
    Returned as a dictionary.

    Note: Here we should use xHermes, but the grid file does not
    include the Y guard cell values. xHermes/xBOUT therefore fails
    to load the dataset due to array size mismatches if Y boundaries
    are requested.

    Since extrapolate_y = false in this case, the zShift angle is
    constant into the boundary. Hence no mapping to field-aligned
    coordinates is needed to reconstruct sheath entrance values.

    NOTE: The code below is specific to this case!
    """

    # Many diagnostics are mapped in flux space to R at midplane
    with DataFile(gridfilepath) as grid:
        psixy = grid["psixy"]
        Rxy = grid["Rxy"]
        ixsep = grid["ixseps1"]

    psi_mid = psixy[:, ymid]
    R_mid = Rxy[:, ymid]
    Rsep = 0.5 * (R_mid[ixsep - 1] + R_mid[ixsep])

    from scipy.interpolate import interp1d

    R_u = interp1d(
        psi_mid, R_mid - Rsep, fill_value="extrapolate"
    )  # Linearly interpolate R as a function of psi

    Nnorm = collect("Nnorm", path=path)
    Tnorm = collect("Tnorm", path=path)
    wci = collect("Omega_ci", path=path)
    Cs0 = collect("Cs0", path=path)
    time = collect("t_array", path=path) / wci  # Seconds
    run_id = collect("run_id", path=path)
    AA = 2  # Atomic species

    # Dictionary to be populated by the add_var function
    result = {}

    def add_var(name, units, data_txyz):
        assert len(data_txyz.shape) == 4
        assert data_txyz.shape[1] == psixy.shape[0]  # X axis
        assert data_txyz.shape[2] == psixy.shape[1] + 4  # Y axis.
        # Data includes 4 Y guards that are not in the grid file

        result[name] = {}

        def add_location(location, name, Rx, data_txz):
            assert len(Rx.shape) == 1
            assert len(data_txz.shape) == 3
            assert len(Rx) == data_txz.shape[1]

            result[name][location] = {
                "name": name,
                "units": units,
                "nt": data_txz.shape[0],
                "tmin": time[0],
                "tmax": time[-1],
                "duration": time[-1] - time[0],
                "nz": data_txz.shape[-1],
                "run_ids": [run_id],
                "paths": [path],
                "gridfilepath": gridfilepath,
                "mean": np.mean(data_txz, axis=(0, -1)),  # Average in time and Z
                "meansq": np.mean(data_txz**2, axis=(0, -1)),
                "std": np.std(data_txz, axis=(0, -1)),
                "Ru": Rx * 100.0,  # Major radius in cm
                "Ru_units": "cm",
            }

        # Outboard midplane. data includes 2 Y guard cells
        add_location("omp", name, R_mid - Rsep, data_txyz[:, :, ymid + 2, :])

        # Interpolate at high field side, mapping to outboard midplane R
        add_location(
            "hfs",
            name,
            R_u(psixy[:, 0]),
            0.5 * (data_txyz[:, :, 1, :] + data_txyz[:, :, 2, :]),
        )

        # Interpolate at low field side
        add_location(
            "lfs",
            name,
            R_u(psixy[:, -1]),
            0.5 * (data_txyz[:, :, -2, :] + data_txyz[:, :, -3, :]),
        )

    Ne = collect("Ne", path=path, yguards=True)
    add_var("Ne", "1/m^3", Ne * Nnorm)

    Ne_floor = np.clip(Ne, 1e-5, None)

    Pe = collect("Pe", path=path, yguards=True)
    add_var("Pe", "Pa", Pe * Nnorm * Tnorm * qe)

    Te = Pe / Ne_floor
    add_var("Te", "eV", Te * Tnorm)

    Pi = collect("Pi", path=path, yguards=True)
    add_var("Pi", "Pa", Pi * Nnorm * Tnorm * qe)

    Ti = Pi / Ne_floor
    add_var("Ti", "eV", Ti * Tnorm)

    phi = collect("phi", path=path, yguards=True)
    add_var("phi", "V", phi * Tnorm)

    Vfl = phi - 3 * Te
    add_var("Vfl", "V", Vfl * Tnorm)

    # Note: NVi and NVe contain mass factors
    NVi = collect("NVi", path=path, yguards=True)
    add_var("Vi", "m/s", NVi / (AA * Ne_floor) * Nnorm * Cs0)

    NVe = collect("NVe", path=path, yguards=True)
    Jpar = NVi / AA - NVe / AA_me
    add_var("Jpar", "A/m^2", Jpar * qe * Nnorm * Cs0)

    return result


def combine_data(dataset1, dataset2):
    """
    Combine two datasets together. Returns a new dataset.
    """
    result = {}
    for name in dataset1.keys():
        result[name] = {}
        for location in dataset1[name].keys():
            data1 = dataset1[name][location]
            data2 = dataset2[name][location]

            # Check that the datasets are compatible
            assert data1["name"] == data2["name"]
            assert data1["units"] == data2["units"]
            assert data1["nz"] == data2["nz"]
            assert data1["gridfilepath"] == data2["gridfilepath"]
            assert data1["Ru_units"] == data2["Ru_units"]

            # Weighting factors for combining means
            nt_tot = data1["nt"] + data2["nt"]
            w1 = data1["nt"] / nt_tot
            w2 = data2["nt"] / nt_tot
            mean = data1["mean"] * w1 + data2["mean"] * w2
            meansq = data1["meansq"] * w1 + data2["meansq"] * w2
            # Re-calculate standard deviation
            std = np.sqrt(meansq - mean**2)

            result[name][location] = {
                "name": data1["name"],
                "units": data1["units"],
                "nt": nt_tot,
                "tmin": min([data1["tmin"], data2["tmin"]]),
                "tmax": max([data1["tmax"], data2["tmax"]]),
                "duration": data1["duration"] + data2["duration"],
                "nz": data1["nz"],
                "run_ids": data1["run_ids"] + data2["run_ids"],
                "paths": data1["paths"] + data2["paths"],
                "gridfilepath": data1["gridfilepath"],
                "mean": mean,
                "meansq": meansq,
                "std": std,
                "Ru": data1["Ru"],
                "Ru_units": data1["Ru_units"],
            }
    return result


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Gather Hermes-3 data from one or more directories"
    )

    parser.add_argument(
        "paths", type=str, nargs="+", help="Paths containing datasets to concatenate"
    )

    parser.add_argument(
        "-o", "--output", default="data.pickle", help="Pickle file to write to"
    )

    parser.add_argument(
        "-g", "--grid", default="bout.grd.nc", type=str, help="Grid file"
    )

    args = parser.parse_args()

    print(f"Got {len(args.paths)} files: {args.paths}")
    print(f"Outputting to '{args.output}'")

    data = extract_data(args.paths[0], args.grid, ymid=18)
    for path in args.paths[1:]:
        data2 = extract_data(path, args.grid, ymid=18)
        data = combine_data(data, data2)

    with open(args.output, "wb") as f:
        pickle.dump(data, f)
