"""Evaluate the arithmetic BOUT++ expression subset on saved slab coordinates."""

import re

import numpy as np
import xarray as xr

from .dataset import AuditError


def evaluate(run, key, time):
    """Evaluate a resolved expression; time is BOUT++ normalised time.

    Only arithmetic expressions and explicitly provided NumPy functions are
    accepted. Named option expressions are expanded by boutdata first.
    """
    options, ds = run.data.options, run.data
    section = options
    for name in key.split(":")[:-1]:
        section = section[name]
    expression = section._substitute_expressions(key.split(":")[-1]).replace("^", "**")
    nx = ds.sizes["x"] - 2 * run.mxg
    ny = ds.sizes["y"] - 2 * run.myg
    symmetric_x = (
        options.get_bool("mesh:symmetricGlobalX")
        if "mesh:symmetricGlobalX" in options
        else False
    )
    symmetric_y = (
        options.get_bool("mesh:symmetricGlobalY")
        if "mesh:symmetricGlobalY" in options
        else False
    )
    # FieldFactory's default includes boundary points; symmetric meshes use
    # physical cell centres. MMS fixtures request symmetric coordinates.
    x = (np.arange(ds.sizes["x"]) - run.mxg + (0.5 if symmetric_x else 0)) / (
        nx if symmetric_x else nx - 1
    )
    y = (
        2
        * np.pi
        * (np.arange(ds.sizes["y"]) - run.myg + (0.5 if symmetric_y else 0))
        / ny
    )
    z = 2 * np.pi * np.arange(ds.sizes["z"]) / ds.sizes["z"]
    names = {
        name: getattr(np, name)
        for name in ("sin", "cos", "tan", "tanh", "exp", "log", "sqrt", "abs", "where")
    }
    names.update(
        x=x[:, None, None],
        y=y[None, :, None],
        z=z[None, None, :],
        t=float(time),
        pi=np.pi,
    )
    # Exclude object access, strings and Python statements from the expression.
    if re.search(r"[\[\]{};'\"=]|\b__", expression) or re.search(
        r"[A-Za-z_]\w*\s*\.", expression
    ):
        raise AuditError(f"Unsupported expression: {key}")
    try:
        value = eval(expression, {"__builtins__": {}}, names)
        value = np.broadcast_to(value, (ds.sizes["x"], ds.sizes["y"], ds.sizes["z"]))
    except (NameError, ValueError, SyntaxError, TypeError) as error:
        raise AuditError(f"Cannot evaluate {key}: {error}") from None
    return xr.DataArray(
        value, dims=("x", "y", "z"), coords={d: ds[d] for d in ("x", "y", "z")}
    )
