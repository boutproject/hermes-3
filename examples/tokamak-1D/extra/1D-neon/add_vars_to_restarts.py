input_path = "."
input_options = "BOUT.inp"

import os

import numpy as np
from boutdata.data import BoutOptionsFile
from boututils.datafile import DataFile

opts = BoutOptionsFile(input_options)

with DataFile(
    os.path.join(input_path, "BOUT.restart.0.nc"), write=True, create=False
) as r0:
    var = np.zeros(r0.size("Pe"))  # Variables should be the same shape as Pe

    for varname in opts.sections():
        if ("function" in opts[varname]) and (varname not in r0):
            # This should be an evolving variable,
            # but is not in the restart file
            try:
                value = opts[varname].evaluate_scalar("function")
            except ValueError, TypeError, KeyError:
                print(
                    f"Failed to evaluate value for {varname}: {opts[varname]['function']}"
                )
                continue
            print(f"Writing field {varname} to restart file with value {value}")
            r0.write(varname, var + value)

            if varname[0] == "N" or varname[0] == "P":
                # Also write the logarithm
                print(
                    "Writing field {} to restart file with value {}".format(
                        "log" + varname, np.log(value)
                    )
                )
                r0.write("log" + varname, var + np.log(value))

        elif (
            ((varname[0] == "N") and (varname[1] != "V")) or (varname[0] == "P")
        ) and ("log" + varname not in r0):
            # Density or pressure but no logarithm
            value = r0[varname]
            minvalue = np.amin(value)
            if minvalue < 1e-10:
                print(f"Minimum value of {varname} is {minvalue}. Flooring")
                value = np.clip(value, 1e-10, None)
                r0.write(varname, value)
            print(
                "Writing field {} with logarithm of field {}".format(
                    "log" + varname, varname
                )
            )
            r0.write("log" + varname, np.log(value))
