from pathlib import Path
import tcvx21
from tcvx21.record_c.record_writer_m import RecordWriter


def write_x21_dataset(hermes_data: dict, output_file: Path):
    """
    Write data to NetCDF, in the format of the TCV-X21 datasets
    """

    result = {
        "LFS-LP": {
            "name": "Low-field-side target Langmuir probes",
            "hermes_location": "lfs",
            "observables": {
                "density": {
                    "name": "Plasma density",
                    "hermes_name": "Ne",
                    "experimental_hierarchy": 2,
                    "dimensionality": 1,
                    "simulation_hierarchy": 1,
                },
                "electron_temp": {
                    "name": "Electron temperature",
                    "hermes_name": "Te",
                    "experimental_hierarchy": 2,
                    "dimensionality": 1,
                    "simulation_hierarchy": 1,
                },
                "ion_temp": {
                    "name": "Ion temperature",
                    "hermes_name": "Ti",
                    "experimental_hierarchy": 2,
                    "dimensionality": 1,
                    "simulation_hierarchy": 1,
                },
                "potential": {
                    "name": "Plasma potential",
                    "hermes_name": "phi",
                    "experimental_hierarchy": 2,
                    "dimensionality": 1,
                    "simulation_hierarchy": 1,
                },
                "current": {
                    "name": "Parallel current",
                    "hermes_name": "Jpar",
                    "experimental_hierarchy": 1,
                    "dimensionality": 1,
                    "simulation_hierarchy": 1,
                },
                "vfloat": {
                    "name": "Floating potential",
                    "hermes_name": "Vfl",
                    "experimental_hierarchy": 1,
                    "dimensionality": 1,
                    "simulation_hierarchy": 2,
                },
            },
        },
        "HFS-LP": {
            "name": "High-field-side target Langmuir probes",
            "hermes_location": "hfs",
            "observables": {
                "density": {
                    "name": "Plasma density",
                    "hermes_name": "Ne",
                    "experimental_hierarchy": 2,
                    "dimensionality": 1,
                    "simulation_hierarchy": 1,
                },
                "electron_temp": {
                    "name": "Electron temperature",
                    "hermes_name": "Te",
                    "experimental_hierarchy": 2,
                    "dimensionality": 1,
                    "simulation_hierarchy": 1,
                },
                "ion_temp": {
                    "name": "Ion temperature",
                    "hermes_name": "Ti",
                    "experimental_hierarchy": 2,
                    "dimensionality": 1,
                    "simulation_hierarchy": 1,
                },
                "potential": {
                    "name": "Plasma potential",
                    "hermes_name": "phi",
                    "experimental_hierarchy": 2,
                    "dimensionality": 1,
                    "simulation_hierarchy": 1,
                },
                "current": {
                    "name": "Parallel current",
                    "hermes_name": "Jpar",
                    "experimental_hierarchy": 1,
                    "dimensionality": 1,
                    "simulation_hierarchy": 1,
                },
                "vfloat": {
                    "name": "Floating potential",
                    "hermes_name": "Vfl",
                    "experimental_hierarchy": 1,
                    "dimensionality": 1,
                    "simulation_hierarchy": 2,
                },
            },
        },
        "FHRP": {
            "name": "Outboard midplane reciprocating probe",
            "hermes_location": "omp",
            "observables": {
                "density": {
                    "name": "Plasma density",
                    "hermes_name": "Ne",
                    "experimental_hierarchy": 2,
                    "dimensionality": 1,
                    "simulation_hierarchy": 1,
                },
                "electron_temp": {
                    "name": "Electron temperature",
                    "hermes_name": "Te",
                    "experimental_hierarchy": 2,
                    "dimensionality": 1,
                    "simulation_hierarchy": 1,
                },
                "ion_temp": {
                    "name": "Ion temperature",
                    "hermes_name": "Ti",
                    "experimental_hierarchy": -1,
                    "dimensionality": 1,
                    "simulation_hierarchy": 1,
                },
                "potential": {
                    "name": "Plasma potential",
                    "hermes_name": "phi",
                    "experimental_hierarchy": 2,
                    "dimensionality": 1,
                    "simulation_hierarchy": 1,
                },
                "vfloat": {
                    "name": "Floating potential",
                    "hermes_name": "Vfl",
                    "experimental_hierarchy": 1,
                    "dimensionality": 1,
                    "simulation_hierarchy": 2,
                },
            },
        },
    }

    # Add Hermes-3 data
    for dname, diagnostic in result.items():
        observables = diagnostic["observables"]
        location = diagnostic["hermes_location"]
        for oname, observable in observables.items():
            hermes_obs = hermes_data[observable["hermes_name"]][location]
            observable["units"] = hermes_obs["units"]
            observable["values"] = hermes_obs["mean"]
            observable["errors"] = hermes_obs["std"]
            observable["Ru"] = hermes_obs["Ru"]
            observable["Ru_units"] = hermes_obs["Ru_units"]

    additional_attributes = {}

    writer = RecordWriter(
        file_path=output_file,
        descriptor="Hermes-3",
        description="Hermes-3 simulation dataset",
        allow_overwrite=True,
    )
    writer.write_data_dict(result, additional_attributes)


if __name__ == "__main__":
    import argparse
    import pickle

    parser = argparse.ArgumentParser(
        description="Convert Hermes-3 pickle file into TCV-X21 NetCDF file"
    )

    parser.add_argument(
        "pickle_file_path", type=str, help="Pickle file containing Hermes-3 data"
    )

    parser.add_argument(
        "-o",
        "--output",
        default="hermes_tcvx21_data.nc",
        help="Output NetCDF file in TCV-X21 format",
    )

    args = parser.parse_args()

    with open(args.pickle_file_path, "rb") as f:
        data = pickle.load(f)

    write_x21_dataset(data, Path(args.output))
