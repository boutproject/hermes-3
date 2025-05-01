import json
import re
import os.path


def construct_inpath(lbl):
    return os.path.abspath(f"./src/amjuel_{lbl}.cxx")


def construct_outpath(lbl):
    return os.path.abspath(f"./json_database/{lbl}.json")


def extract_from_file(fpath, start_line, end_line):
    # Read src
    with open(fpath, "r") as fh:
        lines = fh.readlines()

    # Single string with newlines removed
    s = ("".join(lines[start_line - 1 : end_line])).replace("\n", "")

    # Find rows inside {}
    rows = re.findall(r"\{(.*?)\}", s)
    # Split comma-separated vals and convert to floats
    vals = []
    for row in rows:
        row_vals = [float(s) for s in row.split(",")]
        vals.append(row_vals)
    return vals


def write_json(fpath, rate_coeffs, radiation_coeffs, electron_heating, **meta):
    data = dict(
        info=meta,
        electron_heating=electron_heating,
        rate_coeffs=rate_coeffs,
        radiation_coeffs=radiation_coeffs,
    )
    with open(fpath, "w") as fh:
        json.dump(data, fh)


def gen_data(
    inlbl,
    rate_coeffs_line_range,
    radiation_coeffs_line_range,
    electron_heating,
    outlbl=None,
    rates_amjuel_id="?.?.?",
    rates_amjuel_page="?",
    rad_amjuel_id="?.?.?",
    rad_amjuel_page="?",
    reaction_lbl="?",
    rad_notes="",
    rates_notes="",
):
    # Name json files in the same way as the source files by default
    if outlbl is None:
        outlbl = inlbl

    inpath = construct_inpath(inlbl)
    outpath = construct_outpath(outlbl)
    rate_coeffs = extract_from_file(inpath, *rate_coeffs_line_range)
    radiation_coeffs = extract_from_file(inpath, *radiation_coeffs_line_range)

    # Data labels, etc.
    rate_coeffs_ref = f"Amjuel {rates_amjuel_id}"
    if rates_amjuel_page is not None:
        rate_coeffs_ref += f", page {rates_amjuel_page}"

    rad_coeffs_ref = f"Amjuel {rad_amjuel_id}"
    if rad_amjuel_page is not None:
        rad_coeffs_ref += f", page {rad_amjuel_page}"
    meta = dict(
        rate_coeffs_ref=rate_coeffs_ref,
        rad_coeffs_ref=rad_coeffs_ref,
        reaction_lbl=reaction_lbl,
        rad_notes=rad_notes,
        rates_notes=rates_notes,
    )

    # Write json
    write_json(outpath, rate_coeffs, radiation_coeffs, electron_heating, **meta)


# Extract data and generate jsons for H, He ionization and recombination
gen_data(
    "hyd_ionisation",
    (7, 33),
    (39, 65),
    0.0,
    rates_amjuel_id="H.4 2.1.5",
    rates_amjuel_page="135",
    rad_amjuel_id="H.10 2.1.5",
    rad_amjuel_page="280",
    reaction_lbl="h + e -> h+ + 2e",
    rad_notes="Electron energy loss weighted rate coefficient. Data: Sawada/Fujimoto",
    rates_notes="Effective hydrogenic ionization rate Data: K. Sawada/T. Fujimoto.",
)

gen_data(
    "hyd_recombination",
    (7, 33),
    (39, 65),
    13.6,
    rates_amjuel_id="H.4 2.1.8",
    rates_amjuel_page="141",
    rad_amjuel_id="H.10 2.1.8",
    rad_amjuel_page="284",
    reaction_lbl="h+ + e -> h",
    rad_notes="effective electron cooling rate due to rad.+three-b. recombination potential energy loss 13.6*(effrec.rate) still needs to be subtracted (may render the loss negative, i.e., turn it into a gain) Hence: the quantitity given here happens to be the radiation loss.",
    rates_notes="Effective hydrogenic recombination rate Data: K. Sawada, T.Fujimoto, radiative + three-body contribution",
)

gen_data(
    "helium",
    (9, 35),
    (41, 67),
    0.0,
    outlbl="he_ionisation01",
    rates_amjuel_id="H.4 2.3.9a",
    rates_amjuel_page="161",
    rad_amjuel_id="H.10 2.3.9a",
    rad_amjuel_page="293",
    reaction_lbl="e + he -> he+ + 2e",
    rad_notes="Eth=24.588 eV effective electron cooling rate due to ionization of Helium atoms. Fujimoto Formulation II (only ground level transported, no meta-stables kept explicit)",
    rates_notes="Helium multi-step model, here ionization, Eth=24.56 eV. Fujimoto Formulation II, meta-stable unresolved, (only ground level transported, no metastables kept explicit).",
)

gen_data(
    "helium",
    (88, 114),
    (120, 146),
    24.586,
    outlbl="he_recombination10",
    rates_amjuel_id="H.4 2.3.13a",
    rates_amjuel_page="181",
    rad_amjuel_id="H.10 2.3.13a",
    rad_amjuel_page="295",
    reaction_lbl="e + he+ -> he",
    rad_notes="Helium multi-step model, here recombination: radiative + threebody + dielectronic. Fujimoto Formulation II (only ground level transported, no meta-stables kept explicit). The quantity given here happens to be the radiation loss. The loss of potential energy still needs to be subtracted to make this a total electron energy loss (or gain) rate.",
    rates_notes="Helium multi-step model, here recombination: radiative + threebody + dielectronic. Fujimoto Formulation II (only ground level transported, no metastables kept explicit).",
)
