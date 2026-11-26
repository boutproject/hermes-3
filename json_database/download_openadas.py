#!/usr/bin/env python3
"""Download OpenADAS ADF11 files and convert them to the Hermes-3 JSON format.

Some OpenADAS data used by Hermes-3 (currently argon) are not included in the
repository. This script downloads the ADF11 files from https://open.adas.ac.uk
and writes one JSON file per ADF11 file, in the format read by OpenADAS in
src/adas_reaction.cxx.

Examples (from the Hermes-3 source directory):

    python3 json_database/download_openadas.py ar
    python3 json_database/download_openadas.py scd89_ar plt89_ar
    python3 json_database/download_openadas.py --output-dir my_json_db ar

Without arguments, all files in HERMES_FILES that are missing from the output
directory are created. Existing files are only replaced with --overwrite.

The conversion follows build_json.py of OpenADAS_to_JSON by Thomas Body,
https://github.com/TBody/OpenADAS_to_JSON (Copyright (c) 2017 Tom Body,
MIT License), which produced the OpenADAS JSON files in json_database:
log_density is converted from log10(cm^-3) to log10(m^-3) and log_coeff from
cm^3 to m^3 (i.e. -6 in log10); for the ecd class (ionisation potentials in
eV) the first block (Z1 = 0) is dropped and log10 of the values is stored.
Only the Python standard library is used.
"""

import argparse
import http.client
import json
import math
import re
import sys
import urllib.request
from pathlib import Path

URL = "https://open.adas.ac.uk/download/adf11/{prefix}/{stem}.dat"

# ADF11 files used by Hermes-3, per element: class and two-digit year
HERMES_FILES = {
    "c": "acd96 ccd96 ecd96 plt96 prb96 prc96 scd96",
    "li": "acd96 ccd89 ecd89 ecd96 plt96 prb96 prc89 scd96",
    "n": "acd96 ecd96 plt96 prb96 scd96",
    "ne": "acd96 ccd89 ecd96 plt96 prb96 scd96",
    "ar": "acd89 plt89 prb89 scd89",
}

# ADF11 classes that can be converted (as in OpenADAS_to_JSON)
CLASSES = ("acd", "scd", "ccd", "prb", "plt", "prc", "ecd")

HELP = (
    "JSON file corresponding to an OpenADAS data file\n"
    "Created by hermes-3/json_database/download_openadas.py\n"
    "Format of https://github.com/TBody/OpenADAS_to_JSON"
)


def read_values(line):
    """Return the numbers in a line of 10-character Fortran fields."""
    fields = (line[i : i + 10] for i in range(0, len(line), 10))
    return [float(field) for field in fields if field.strip()]


def read_adf11(text):
    """Read an unresolved ADF11 file.

    Returns the nuclear charge, log10(density [cm^-3]), log10(temperature [eV])
    and the data blocks, each a list over temperature of lists over density.
    """
    lines = text.splitlines()
    try:
        charge, n_dens, n_temp = (int(value) for value in lines[0].split()[:3])
    except (IndexError, ValueError):
        raise ValueError("not an ADF11 file") from None
    pos = 2  # Skip the header and the separator line
    grid = []
    while len(grid) < n_dens + n_temp:
        grid += read_values(lines[pos])
        pos += 1
    blocks = []
    # Data blocks start with a line containing "Z1="; comments start with "C"
    while pos < len(lines) and not lines[pos].startswith("C"):
        if "Z1=" not in lines[pos]:
            pos += 1
            continue
        pos += 1
        values = []
        while len(values) < n_dens * n_temp:
            values += read_values(lines[pos])
            pos += 1
        blocks.append([values[i * n_dens : (i + 1) * n_dens] for i in range(n_temp)])
    if not blocks:
        raise ValueError("no data blocks found")
    return charge, grid[:n_dens], grid[n_dens : n_dens + n_temp], blocks


def convert(text, stem, name):
    """Convert the contents of an ADF11 file to the JSON text of OpenADAS_to_JSON."""
    adf11_class, element = stem[:3], stem.split("_")[1]
    charge, log_density, log_temperature, blocks = read_adf11(text)
    if adf11_class == "ecd":
        log_coeff = [
            [[math.log10(value) for value in row] for row in block]
            for block in blocks[1:]
        ]
    else:
        log_coeff = [
            [[value - 6 for value in row] for row in block] for block in blocks
        ]
    data = {
        "charge": charge,
        "class": adf11_class,
        "element": element,
        "help": HELP,
        "log_coeff": log_coeff,
        "log_density": [value + 6 for value in log_density],
        "log_temperature": log_temperature,
        "name": name,
        "number_of_charge_states": len(blocks),
        "numpy_ndarrays": ["log_density", "log_temperature", "log_coeff"],
    }
    return json.dumps(data, sort_keys=True, indent=4)


def file_names(names):
    """Expand element symbols to file names, check file names."""
    if not names:
        names = list(HERMES_FILES)
    stems = []
    for name in names:
        if name in HERMES_FILES:
            stems += [f"{prefix}_{name}" for prefix in HERMES_FILES[name].split()]
        elif re.fullmatch(r"[a-z]{3}[0-9]{2}_[a-z]{1,2}", name) and name[:3] in CLASSES:
            stems.append(name)
        else:
            raise SystemExit(
                f"Unknown name '{name}': expected an element in {list(HERMES_FILES)}"
                f" or an ADF11 file name such as scd89_ar (classes {CLASSES})"
            )
    return stems


def download(url):
    """Download an ADF11 file and return its contents."""
    with urllib.request.urlopen(url, timeout=60) as response:
        content_type = response.headers.get("Content-Type", "")
        data = response.read()
    if not content_type.startswith("data/adf11"):
        raise ValueError(f"not an ADF11 file (Content-Type '{content_type}')")
    return data


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "names",
        nargs="*",
        help="elements (e.g. ar) or ADF11 file names (e.g. scd89_ar); default: all",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="directory for the JSON files (default: %(default)s)",
    )
    parser.add_argument(
        "--dat-dir",
        type=Path,
        help="keep the downloaded .dat files here, and use files already present",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="replace existing JSON files"
    )
    args = parser.parse_args(argv)

    failed = []
    for stem in file_names(args.names):
        output = args.output_dir / f"{stem}.json"
        if output.exists() and not args.overwrite:
            print(f"{output} exists, skipped (use --overwrite to replace it)")
            continue
        url = URL.format(prefix=stem.split("_")[0], stem=stem)
        try:
            dat_file = args.dat_dir / f"{stem}.dat" if args.dat_dir else None
            if dat_file and dat_file.exists():
                data = dat_file.read_bytes()
            else:
                data = download(url)
                if dat_file:
                    dat_file.parent.mkdir(parents=True, exist_ok=True)
                    dat_file.write_bytes(data)
            json_text = convert(data.decode("latin-1"), stem, url)
        # http.client.HTTPException: e.g. IncompleteRead if a download is interrupted
        except (OSError, ValueError, IndexError, http.client.HTTPException) as error:
            print(f"{url}: {error}", file=sys.stderr)
            failed.append(stem)
            continue
        output.parent.mkdir(parents=True, exist_ok=True)
        with open(output, "w", encoding="ascii") as fp:
            fp.write(json_text)
        print(f"{url} -> {output}")
    if failed:
        raise SystemExit(f"Failed: {' '.join(failed)}")


if __name__ == "__main__":
    main()
