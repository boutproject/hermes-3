#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import os
import re
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path


def setup_plot_environment() -> None:
    if not any(arg in {"--plot", "--show", "--anim"} for arg in sys.argv[1:]):
        os.environ.setdefault("MPLBACKEND", "Agg")
    mplconfigdir = Path(os.environ.get("TMPDIR", "/tmp")) / "hermes_ib_mplconfig"
    mplconfigdir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mplconfigdir))


setup_plot_environment()

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from boutdata import collect
from boututils.datafile import DataFile
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.widgets import Button, Slider
from scipy import special
from scipy.interpolate import RectBivariateSpline


HERE = Path(__file__).resolve().parent
DEFAULT_OMEGA_CI = 95788333.030660808
BC_RE = re.compile(
    r"^(?P<kind>dirichlet|neumann|value)\((?P<value>[-+0-9.eE]+)\)$",
    re.IGNORECASE,
)
SAVE_ANIM_DEFAULT = "__default_animation_path__"


@dataclass
class GhostError:
    time_index: int
    time_seconds: float
    time_norm: float
    n_ghost: int
    l2: float
    linf: float


@dataclass
class RunResult:
    bc_name: str
    resolution: int
    case_dir: Path
    grid_path: Path
    h: float
    density: np.ndarray
    truth: np.ndarray
    error: np.ndarray
    mask: np.ndarray
    ghost_mask: np.ndarray
    mxg: int
    R: np.ndarray
    Z: np.ndarray
    X: np.ndarray
    Zn: np.ndarray
    l2: float
    linf: float
    function_expr: str
    lambda_value: float | None
    time_model: str
    time_index: int
    time_seconds: float
    time_norm: float
    ghost_t0: GhostError
    ghost_selected: GhostError
    spline_density: RectBivariateSpline
    spline_truth: RectBivariateSpline


@dataclass
class InterpRun:
    bc_name: str
    resolution: int
    case_dir: Path
    grid_path: Path
    h: float
    density: np.ndarray
    mask: np.ndarray
    ghost_mask: np.ndarray
    mxg: int
    R: np.ndarray
    Z: np.ndarray
    t_array: np.ndarray
    t_norm: np.ndarray
    time_index: int
    time_seconds: float
    time_norm_selected: float


@dataclass
class InterpPair:
    bc_name: str
    coarse: InterpRun
    fine: InterpRun
    h: float
    mask: np.ndarray
    point_count: int
    l2: float
    linf: float
    l2_series: np.ndarray
    linf_series: np.ndarray
    error_selected: np.ndarray

    @property
    def label(self) -> str:
        return f"{self.coarse.resolution} -> {self.fine.resolution}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze IB MMS density convergence for dirichlet/neumann runs."
    )
    parser.add_argument(
        "directory",
        nargs="?",
        default=None,
        help="Case directory or case name under examples/ib/mms, e.g. diamond or circle.",
    )
    parser.add_argument(
        "--d",
        "--case-dir",
        "--case",
        dest="directory_option",
        default=None,
        help="Case directory or case name under examples/ib/mms, matching grid_gen.py.",
    )
    parser.add_argument(
        "--dir",
        dest="plot_dirichlet",
        action="store_true",
        help="Use Dirichlet cases for the single-family density, error, setup, and animation plots.",
    )
    parser.add_argument(
        "--tind",
        dest="time_index",
        type=int,
        default=None,
        help=(
            "Time index to analyze and plot. Defaults to the halfway timestep; "
            "negative indices follow Python indexing."
        ),
    )
    parser.add_argument(
        "--time-index",
        dest="time_index",
        type=int,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Density/convergence PDF path. Defaults to ./out/<case>MMS.pdf.",
    )
    parser.add_argument(
        "--ghost-output",
        type=Path,
        default=None,
        help="Ghost convergence PDF path. Defaults to ./out/<case>GhostMMS.pdf.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("out"),
        help="Directory for default PDF outputs. Defaults to ./out.",
    )
    parser.add_argument(
        "--plot",
        "--show",
        dest="plot",
        action="store_true",
        help="Show plots interactively after saving.",
    )
    parser.add_argument(
        "--anim",
        action="store_true",
        help="Show an animated highest-resolution density/error diagnostic.",
    )
    parser.add_argument(
        "--save-anim",
        nargs="?",
        const=SAVE_ANIM_DEFAULT,
        default=None,
        metavar="PATH",
        help=(
            "Save the highest-resolution density/error animation. If PATH is "
            "omitted, writes ./out/<case>Anim.gif or ./out/<case>InterpAnim.gif."
        ),
    )
    parser.add_argument(
        "--interp",
        action="store_true",
        help=(
            "Compare against the next higher-resolution run instead of the analytic "
            "MMS solution. This is intended for non-analytic cases such as d3d."
        ),
    )
    parser.add_argument(
        "--interp-order",
        type=int,
        default=1,
        choices=(1, 2, 3),
        help="Spline order for --interp successive comparisons. Defaults to 1.",
    )
    parser.add_argument(
        "--log",
        dest="log_error",
        action="store_true",
        help="Plot spatial error panels and error animations as log10(error).",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=200,
        help="DPI used for raster outputs.",
    )
    parser.add_argument(
        "--plot-points",
        type=int,
        default=300,
        help=(
            "Deprecated compatibility option. Analytic density/error panels "
            "are plotted on the raw masked grid."
        ),
    )
    return parser.parse_args()


def resolve_case_dir(directory: str) -> Path:
    candidate = Path(directory).expanduser()
    if candidate.exists():
        return candidate.resolve()

    base = Path(__file__).resolve().parent
    candidate = base / directory
    if candidate.exists():
        return candidate.resolve()

    raise FileNotFoundError(f"Could not find MMS case directory: {directory}")


def parse_bout_inp(path: Path) -> dict[str, dict[str, str]]:
    sections: dict[str, dict[str, str]] = {}
    current: str | None = None
    for raw in path.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            current = line[1:-1].strip()
            sections.setdefault(current, {})
            continue
        if "=" in line and current is not None:
            key, value = [item.strip() for item in line.split("=", 1)]
            sections[current][key] = value.strip().strip('"')
    return sections


def parse_boundary_condition(expr: str) -> tuple[str, float]:
    match = BC_RE.fullmatch(expr.replace(" ", ""))
    if match is None:
        raise ValueError(f"Unsupported boundary condition: {expr}")
    kind = match.group("kind").lower()
    if kind == "value":
        kind = "dirichlet"
    return kind, float(match.group("value"))


def canonical_bc_name(name: str) -> str:
    aliases = {
        "dir": "dirichlet",
        "dirichlet": "dirichlet",
        "neu": "neumann",
        "neumann": "neumann",
    }
    return aliases.get(name.lower(), name.lower())


def sympy_expression(expr: str) -> sp.Expr:
    x, z, zn = sp.symbols("x z zn")

    def _besselj(order, arg):
        return sp.besselj(order, arg)

    local_dict = {
        "x": x,
        "z": z,
        "zn": zn,
        "pi": sp.pi,
        "sin": sp.sin,
        "cos": sp.cos,
        "sqrt": sp.sqrt,
        "exp": sp.exp,
        "j": _besselj,
        "besselj": sp.besselj,
    }
    parsed = sp.sympify(expr.replace("^", "**"), locals=local_dict)
    return sp.simplify(parsed.subs(zn, z))


def split_constant_spatial(expr: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    x, z = sp.symbols("x z")
    expr = sp.expand(expr)
    constant = sp.Integer(0)
    for term in sp.Add.make_args(expr):
        if not (term.free_symbols & {x, z}):
            constant += term
    spatial = sp.simplify(expr - constant)
    return sp.simplify(constant), sp.simplify(spatial)


def lambdify_expr(expr: sp.Expr):
    x, z = sp.symbols("x z")
    modules = [{"besselj": special.jv}, "numpy"]
    return sp.lambdify((x, z), expr, modules=modules)


def lambdify_time_expr(expr: sp.Expr):
    x, z, t = sp.symbols("x z t")
    modules = [{"besselj": special.jv}, "numpy"]
    return sp.lambdify((x, z, t), expr, modules=modules)


def infer_eigenvalue(
    spatial: sp.Expr,
    X: np.ndarray,
    Zn: np.ndarray,
    mask: np.ndarray,
    x_laplacian_scale: float,
    z_laplacian_scale: float,
) -> float | None:
    x, z = sp.symbols("x z")
    if spatial == 0:
        return 0.0

    laplacian = (
        x_laplacian_scale * sp.diff(spatial, x, 2)
        + z_laplacian_scale * sp.diff(spatial, z, 2)
    )
    ratio = sp.simplify(-laplacian / spatial)
    if not (ratio.free_symbols & {x, z}):
        value = complex(sp.N(ratio))
        if abs(value.imag) < 1e-10:
            return float(value.real)

    try:
        ratio_fn = lambdify_expr(ratio)
        values = np.asarray(ratio_fn(X, Zn), dtype=float)
    except Exception:
        return None

    spatial_fn = lambdify_expr(spatial)
    spatial_values = np.asarray(spatial_fn(X, Zn), dtype=float)
    valid = (
        mask
        & np.isfinite(values)
        & np.isfinite(spatial_values)
        & (np.abs(spatial_values) > 1.0e-10)
    )
    if np.count_nonzero(valid) < 4:
        return None

    sample = values[valid]
    mean = float(np.nanmean(sample))
    spread = float(np.nanmax(np.abs(sample - mean)))
    scale = max(1.0, abs(mean))
    if spread / scale < 1.0e-7:
        return mean

    return None


def is_zero_expr(expr: sp.Expr) -> bool:
    expr = sp.simplify(expr)
    return expr == 0 or bool(expr.equals(0))


def finite_heat_series(
    expr: sp.Expr,
    D: float,
    x_laplacian_scale: float,
    z_laplacian_scale: float,
    max_terms: int = 12,
) -> sp.Expr | None:
    x, z, t = sp.symbols("x z t")
    laplacian = sp.simplify(expr)
    terms: list[sp.Expr] = []

    for order in range(max_terms):
        if is_zero_expr(laplacian):
            return sp.simplify(sum(terms, sp.Integer(0)))
        terms.append((D * t) ** order * laplacian / sp.factorial(order))
        laplacian = sp.simplify(
            x_laplacian_scale * sp.diff(laplacian, x, 2)
            + z_laplacian_scale * sp.diff(laplacian, z, 2)
        )

    return None


def build_time_solution(
    expr: sp.Expr,
    D: float,
    X: np.ndarray,
    Zn: np.ndarray,
    mask: np.ndarray,
    x_laplacian_scale: float,
    z_laplacian_scale: float,
):
    x, z, t = sp.symbols("x z t")
    constant_expr, spatial_expr = split_constant_spatial(expr)
    lambda_value = infer_eigenvalue(
        spatial_expr,
        X,
        Zn,
        mask,
        x_laplacian_scale,
        z_laplacian_scale,
    )
    if lambda_value is not None:
        time_expr = constant_expr + spatial_expr * sp.exp(-D * lambda_value * t)
        return lambdify_time_expr(time_expr), lambda_value, "single eigenmode"

    finite_expr = finite_heat_series(expr, D, x_laplacian_scale, z_laplacian_scale)
    if finite_expr is not None:
        return lambdify_time_expr(finite_expr), None, "finite heat series"

    warnings.warn(
        "Could not derive time dependence for the input function; using the "
        "t=0 function at all times for this run."
    )
    static_fn = lambdify_expr(expr)

    def static_solution(R_value, Z_value, _t_value):
        return static_fn(R_value, Z_value)

    return static_solution, None, "static t=0 fallback"


def grid_values(value, shape: tuple[int, ...]) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.shape == ():
        return np.full(shape, float(arr))
    return np.broadcast_to(arr, shape).astype(float, copy=False)


def datafile_scalar(path: Path, name: str, default: float | None = None) -> float | None:
    try:
        with DataFile(str(path)) as datafile:
            value = datafile.read(name)
    except Exception:
        return default
    arr = np.asarray(value)
    if arr.size == 0:
        return default
    return float(arr.reshape(-1)[0])


def first_dump_file(case_dir: Path) -> Path:
    dumps = sorted(case_dir.glob("BOUT.dmp.*.nc"))
    if not dumps:
        raise FileNotFoundError(f"No BOUT.dmp.*.nc files found in {case_dir}")
    return dumps[0]


def read_grid(grid_path: Path) -> dict[str, np.ndarray | float | int]:
    grid: dict[str, np.ndarray | float | int] = {}
    with DataFile(str(grid_path)) as datafile:
        names = ("R", "Z", "in_mask", "vol_frac", "MXG", "nx", "ny", "nz")
        optional = (
            "ghost_id",
            "ghost_pts",
            "image_pts",
            "bndry_pts",
            "vol_frac",
            "cc_mask",
            "bound_id",
            "bound_count",
            "eff_bndry_vertices",
        )
        available = set(datafile.keys())
        for name in names:
            grid[name] = datafile.read(name)
        for name in optional:
            if name in available:
                grid[name] = datafile.read(name)
    return grid


def squeeze_y(arr: np.ndarray) -> np.ndarray:
    arr = np.asarray(arr)
    if arr.ndim == 4:
        return arr[:, :, 0, :]
    if arr.ndim == 3:
        if arr.shape[1] == 1:
            return arr[:, 0, :]
        return arr
    if arr.ndim == 2:
        return arr
    raise ValueError(f"Expected a 2D, 3D, or 4D array, got shape {arr.shape}")


def squeeze_grid(arr: np.ndarray) -> np.ndarray:
    arr = np.asarray(arr)
    if arr.ndim == 3:
        return arr[:, 0, :]
    if arr.ndim == 2:
        return arr
    raise ValueError(f"Expected a 2D or 3D grid array, got shape {arr.shape}")


def selected_time_index(nt: int, time_index: int | None) -> int:
    index = nt // 2 if time_index is None else time_index
    index = index if index >= 0 else nt + index
    if index < 0 or index >= nt:
        raise IndexError(f"time index {time_index} is out of range for nt={nt}")
    return index


def core_domain_mask(shape: tuple[int, int], mxg: int) -> np.ndarray:
    mask = np.ones(shape, dtype=bool)
    if mxg > 0:
        mask[:mxg, :] = False
        mask[-mxg:, :] = False
    return mask


def field_factory_coordinates(shape: tuple[int, int], mxg: int) -> tuple[np.ndarray, np.ndarray]:
    nx_total, nz = shape
    nx_phys = nx_total - 2 * mxg
    if nx_phys <= 0:
        raise ValueError(f"Invalid physical x resolution: nx={nx_total}, MXG={mxg}")
    i = np.arange(nx_total, dtype=float)[:, None]
    k = np.arange(nz, dtype=float)[None, :]
    x = (i - mxg + 0.5) / nx_phys
    zn = (k + 0.5) / nz
    return np.broadcast_to(x, shape), np.broadcast_to(zn, shape)


def error_norms(diff: np.ndarray) -> tuple[float, float]:
    diff = np.asarray(diff, dtype=float)
    finite = np.isfinite(diff)
    if not np.any(finite):
        return math.nan, math.nan
    sample = diff[finite]
    return float(np.sqrt(np.mean(sample**2))), float(np.max(np.abs(sample)))


def centers_to_edges(centers: np.ndarray) -> np.ndarray:
    centers = np.asarray(centers, dtype=float)
    if centers.size < 2:
        dx = 1.0
        return np.array([centers[0] - 0.5 * dx, centers[0] + 0.5 * dx])
    mids = 0.5 * (centers[:-1] + centers[1:])
    edges = np.empty(centers.size + 1, dtype=float)
    edges[1:-1] = mids
    edges[0] = centers[0] - (mids[0] - centers[0])
    edges[-1] = centers[-1] + (centers[-1] - mids[-1])
    return edges


def physical_to_field_coordinates(
    r_values: np.ndarray,
    z_values: np.ndarray,
    R: np.ndarray,
    Z: np.ndarray,
    mxg: int,
) -> tuple[np.ndarray, np.ndarray]:
    x_slice = slice(mxg, -mxg) if mxg > 0 else slice(None)
    r_edges = centers_to_edges(R[x_slice, 0])
    z_edges = centers_to_edges(Z[0, :])
    x = (np.asarray(r_values, dtype=float) - r_edges[0]) / (r_edges[-1] - r_edges[0])
    zn = (np.asarray(z_values, dtype=float) - z_edges[0]) / (z_edges[-1] - z_edges[0])
    return x, zn


def laplacian_coordinate_scales(R: np.ndarray, Z: np.ndarray, mxg: int) -> tuple[float, float]:
    x_slice = slice(mxg, -mxg) if mxg > 0 else slice(None)
    r_edges = centers_to_edges(R[x_slice, 0])
    z_edges = centers_to_edges(Z[0, :])
    r_length = float(r_edges[-1] - r_edges[0])
    z_length = float(z_edges[-1] - z_edges[0])
    if r_length == 0.0 or z_length == 0.0:
        raise ValueError("Cannot infer MMS physical domain length from grid coordinates")
    return 1.0 / (r_length * r_length), 1.0 / (z_length * z_length)


def ghost_error_against_analytic(
    density_t: np.ndarray,
    solution_fn,
    t_array: np.ndarray,
    t_norm: np.ndarray,
    grid: dict[str, np.ndarray | float | int],
    R: np.ndarray,
    Z: np.ndarray,
    time_index: int,
) -> GhostError:
    if "ghost_id" not in grid:
        return GhostError(
            time_index=time_index,
            time_seconds=float(t_array[time_index]),
            time_norm=float(t_norm[time_index]),
            n_ghost=0,
            l2=math.nan,
            linf=math.nan,
        )

    ghost_id = squeeze_grid(np.asarray(grid["ghost_id"]))
    ghost_mask = ghost_id >= 0
    mxg = int(np.asarray(grid["MXG"]).reshape(-1)[0])
    ghost_mask &= core_domain_mask(ghost_mask.shape, mxg)
    if not np.any(ghost_mask):
        return GhostError(
            time_index=time_index,
            time_seconds=float(t_array[time_index]),
            time_norm=float(t_norm[time_index]),
            n_ghost=0,
            l2=math.nan,
            linf=math.nan,
        )

    phi_ghost = np.asarray(density_t, dtype=float)[ghost_mask]
    if "ghost_pts" in grid:
        ghost_pts = np.asarray(grid["ghost_pts"], dtype=float)
        gids = np.rint(ghost_id[ghost_mask]).astype(int)
        valid = (gids >= 0) & (gids < ghost_pts.shape[0])
        phi_ghost = phi_ghost[valid]
        coords = ghost_pts[gids[valid]]
        ghost_R = coords[:, 0]
        ghost_Z = coords[:, 1]
        eval_x, eval_zn = physical_to_field_coordinates(ghost_R, ghost_Z, R, Z, mxg)
    else:
        ghost_R = R[ghost_mask]
        ghost_Z = Z[ghost_mask]
        eval_x, eval_zn = field_factory_coordinates(R.shape, mxg)
        eval_x = eval_x[ghost_mask]
        eval_zn = eval_zn[ghost_mask]

    truth = grid_values(solution_fn(eval_x, eval_zn, t_norm[time_index]), phi_ghost.shape)
    l2, linf = error_norms(phi_ghost - truth)
    return GhostError(
        time_index=time_index,
        time_seconds=float(t_array[time_index]),
        time_norm=float(t_norm[time_index]),
        n_ghost=int(phi_ghost.size),
        l2=l2,
        linf=linf,
    )


def load_case(case_dir: Path, bc_name: str, time_index: int | None) -> RunResult:
    bout = parse_bout_inp(case_dir / "BOUT.inp")
    mesh_file = Path(bout["mesh"]["file"])
    grid_path = mesh_file if mesh_file.is_absolute() else (case_dir / mesh_file).resolve()
    if not grid_path.exists():
        raise FileNotFoundError(f"Grid file for {case_dir} does not exist: {grid_path}")

    function_expr = bout["Nh+"]["function"]
    bc_type, _ = parse_boundary_condition(bout["Nh+"]["bndry_wall"])
    if bc_type != canonical_bc_name(bc_name):
        warnings.warn(f"{case_dir} is under {bc_name}, but BOUT.inp sets {bc_type}.")

    D = float(bout["h+"]["anomalous_D"])
    density = squeeze_y(collect("Nh+", path=str(case_dir), info=False))
    t_array = np.asarray(collect("t_array", path=str(case_dir), info=False), dtype=float)
    nt = density.shape[0]
    t_index = selected_time_index(nt, time_index)

    grid = read_grid(grid_path)
    R = squeeze_grid(np.asarray(grid["R"], dtype=float))
    Z = squeeze_grid(np.asarray(grid["Z"], dtype=float))
    mask = squeeze_grid(np.asarray(grid["in_mask"])).astype(bool)
    MXG = int(np.asarray(grid["MXG"]).reshape(-1)[0])
    if "ghost_id" in grid:
        ghost_mask = squeeze_grid(np.asarray(grid["ghost_id"])) >= 0
    else:
        ghost_mask = np.zeros_like(mask, dtype=bool)
    core_mask = core_domain_mask(mask.shape, MXG)
    mask = mask & core_mask & ~ghost_mask
    ghost_mask = ghost_mask & core_mask
    nx_phys = R.shape[0] - 2 * MXG
    if nx_phys <= 0:
        raise ValueError(f"Invalid physical x resolution for {grid_path}: nx={R.shape[0]}, MXG={MXG}")

    X, Zn = field_factory_coordinates(R.shape, MXG)
    x_laplacian_scale, z_laplacian_scale = laplacian_coordinate_scales(R, Z, MXG)
    expr = sympy_expression(function_expr)
    solution_fn, lambda_value, time_model = build_time_solution(
        expr,
        D,
        X,
        Zn,
        mask,
        x_laplacian_scale,
        z_laplacian_scale,
    )

    omega_ci = datafile_scalar(first_dump_file(case_dir), "Omega_ci", DEFAULT_OMEGA_CI)
    if omega_ci is None or omega_ci == 0.0:
        omega_ci = DEFAULT_OMEGA_CI
    t_norm = t_array / omega_ci
    initial_truth = grid_values(solution_fn(X, Zn, t_norm[0]), R.shape)
    initial_l2, initial_linf = error_norms(np.where(mask, density[0] - initial_truth, np.nan))
    if initial_linf > 1.0e-8:
        warnings.warn(
            f"{case_dir}: t=0 analytic mismatch over plasma cells "
            f"(L2={initial_l2:.3e}, Linf={initial_linf:.3e})."
        )
    truth = grid_values(solution_fn(X, Zn, t_norm[t_index]), R.shape)
    ghost_t0 = ghost_error_against_analytic(
        density[0],
        solution_fn,
        t_array,
        t_norm,
        grid,
        R,
        Z,
        0,
    )
    ghost_selected = ghost_error_against_analytic(
        density[t_index],
        solution_fn,
        t_array,
        t_norm,
        grid,
        R,
        Z,
        t_index,
    )

    density_t = density[t_index]
    diff = density_t - truth
    masked = np.where(mask, diff, np.nan)
    l2 = float(np.sqrt(np.nanmean(masked**2)))
    linf = float(np.nanmax(np.abs(masked)))

    R1d = R[:, 0]
    Z1d = Z[0, :]
    spline_density = RectBivariateSpline(R1d, Z1d, density_t, kx=3, ky=3)
    spline_truth = RectBivariateSpline(R1d, Z1d, truth, kx=3, ky=3)

    resolution_match = re.fullmatch(r"\d+", case_dir.name)
    resolution = int(case_dir.name) if resolution_match else nx_phys
    return RunResult(
        bc_name=bc_name,
        resolution=resolution,
        case_dir=case_dir,
        grid_path=grid_path,
        h=1.0 / float(nx_phys),
        density=density_t,
        truth=truth,
        error=np.abs(diff),
        mask=mask,
        ghost_mask=ghost_mask,
        mxg=MXG,
        R=R,
        Z=Z,
        X=X,
        Zn=Zn,
        l2=l2,
        linf=linf,
        function_expr=function_expr,
        lambda_value=lambda_value,
        time_model=time_model,
        time_index=t_index,
        time_seconds=float(t_array[t_index]),
        time_norm=float(t_norm[t_index]),
        ghost_t0=ghost_t0,
        ghost_selected=ghost_selected,
        spline_density=spline_density,
        spline_truth=spline_truth,
    )


def numeric_case_dirs(base: Path, bc_name: str) -> list[Path]:
    parent = base / bc_name
    if not parent.is_dir():
        return []
    cases = [path for path in parent.iterdir() if path.is_dir() and path.name.isdigit()]
    cases.sort(key=lambda path: int(path.name))
    return cases


def load_family(base: Path, bc_name: str, time_index: int | None) -> list[RunResult]:
    results: list[RunResult] = []
    for case_dir in numeric_case_dirs(base, bc_name):
        if not list(case_dir.glob("BOUT.dmp.*.nc")):
            warnings.warn(f"Skipping {case_dir}: no BOUT.dmp.*.nc files found.")
            continue
        results.append(load_case(case_dir, bc_name, time_index))
    if not results:
        raise FileNotFoundError(f"No completed {bc_name} runs found under {base / bc_name}")
    return results


def fit_order_values(h: np.ndarray, err: np.ndarray) -> float:
    h = np.asarray(h, dtype=float)
    err = np.asarray(err, dtype=float)
    valid = np.isfinite(h) & np.isfinite(err) & (h > 0.0) & (err > 0.0)
    if np.count_nonzero(valid) < 2:
        return math.nan
    return float(np.polyfit(np.log(h[valid]), np.log(err[valid]), 1)[0])


def fit_order(results: list[RunResult], metric: str = "l2") -> float:
    h = np.asarray([run.h for run in results], dtype=float)
    err = np.asarray([getattr(run, metric) for run in results], dtype=float)
    return fit_order_values(h, err)


def bc_display_name(bc_name: str) -> str:
    name = canonical_bc_name(bc_name)
    if name == "dirichlet":
        return "Dirichlet"
    if name == "neumann":
        return "Neumann"
    return bc_name


def selected_plot_family(
    neumann: list[RunResult],
    dirichlet: list[RunResult],
    use_dirichlet: bool,
) -> list[RunResult]:
    return dirichlet if use_dirichlet else neumann


def core_x_slice(run: RunResult) -> slice:
    if run.mxg <= 0:
        return slice(None)
    return slice(run.mxg, -run.mxg)


def load_interp_case(case_dir: Path, bc_name: str, time_index: int | None) -> InterpRun:
    bout = parse_bout_inp(case_dir / "BOUT.inp")
    mesh_file = Path(bout["mesh"]["file"])
    grid_path = mesh_file if mesh_file.is_absolute() else (case_dir / mesh_file).resolve()
    if not grid_path.exists():
        raise FileNotFoundError(f"Grid file for {case_dir} does not exist: {grid_path}")

    density = np.asarray(squeeze_y(collect("Nh+", path=str(case_dir), info=False)), dtype=float)
    t_array = np.asarray(collect("t_array", path=str(case_dir), info=False), dtype=float)
    nt = density.shape[0]
    t_index = selected_time_index(nt, time_index)

    grid = read_grid(grid_path)
    R = squeeze_grid(np.asarray(grid["R"], dtype=float))
    Z = squeeze_grid(np.asarray(grid["Z"], dtype=float))
    mask = squeeze_grid(np.asarray(grid["in_mask"])).astype(bool)
    MXG = int(np.asarray(grid["MXG"]).reshape(-1)[0])
    if "ghost_id" in grid:
        ghost_mask = squeeze_grid(np.asarray(grid["ghost_id"])) >= 0
    else:
        ghost_mask = np.zeros_like(mask, dtype=bool)

    if density.shape[1:] != R.shape:
        raise ValueError(
            f"Density shape {density.shape[1:]} does not match grid shape {R.shape} for {case_dir}"
        )

    core_mask = core_domain_mask(mask.shape, MXG)
    mask = mask & core_mask & ~ghost_mask
    ghost_mask = ghost_mask & core_mask

    nx_phys = R.shape[0] - 2 * MXG
    if nx_phys <= 0:
        raise ValueError(f"Invalid physical x resolution for {grid_path}: nx={R.shape[0]}, MXG={MXG}")

    omega_ci = datafile_scalar(first_dump_file(case_dir), "Omega_ci", DEFAULT_OMEGA_CI)
    if omega_ci is None or omega_ci == 0.0:
        omega_ci = DEFAULT_OMEGA_CI
    t_norm = t_array / omega_ci

    resolution_match = re.fullmatch(r"\d+", case_dir.name)
    resolution = int(case_dir.name) if resolution_match else nx_phys
    return InterpRun(
        bc_name=bc_name,
        resolution=resolution,
        case_dir=case_dir,
        grid_path=grid_path,
        h=1.0 / float(nx_phys),
        density=density,
        mask=mask,
        ghost_mask=ghost_mask,
        mxg=MXG,
        R=R,
        Z=Z,
        t_array=t_array,
        t_norm=t_norm,
        time_index=t_index,
        time_seconds=float(t_array[t_index]),
        time_norm_selected=float(t_norm[t_index]),
    )


def load_interp_family(base: Path, bc_name: str, time_index: int | None) -> list[InterpRun]:
    results: list[InterpRun] = []
    for case_dir in numeric_case_dirs(base, bc_name):
        if not list(case_dir.glob("BOUT.dmp.*.nc")):
            warnings.warn(f"Skipping {case_dir}: no BOUT.dmp.*.nc files found.")
            continue
        results.append(load_interp_case(case_dir, bc_name, time_index))
    if len(results) < 2:
        warnings.warn(
            f"Need at least two completed {bc_name} runs for --interp under {base / bc_name}; "
            f"found {len(results)}."
        )
    return results


def spline_to_target(
    values: np.ndarray,
    r_src: np.ndarray,
    z_src: np.ndarray,
    r_target: np.ndarray,
    z_target: np.ndarray,
    order: int,
) -> np.ndarray:
    kx = min(order, len(r_src) - 1)
    ky = min(order, len(z_src) - 1)
    spline = RectBivariateSpline(
        r_src,
        z_src,
        np.asarray(values, dtype=float),
        kx=kx,
        ky=ky,
    )
    return np.asarray(spline(r_target, z_target), dtype=float)


def interp_overlap_mask(coarse: InterpRun, fine: InterpRun) -> np.ndarray:
    coarse_x_slice = core_x_slice(coarse)
    r_src = coarse.R[coarse_x_slice, 0]
    z_src = coarse.Z[0, :]
    return (
        fine.mask
        & (fine.R >= float(r_src[0]))
        & (fine.R <= float(r_src[-1]))
        & (fine.Z >= float(z_src[0]))
        & (fine.Z <= float(z_src[-1]))
    )


def interp_pair_diff(pair: InterpPair, time_index: int, order: int) -> np.ndarray:
    coarse = pair.coarse
    fine = pair.fine
    coarse_on_fine = spline_to_target(
        coarse.density[time_index],
        coarse.R[:, 0],
        coarse.Z[0, :],
        fine.R[:, 0],
        fine.Z[0, :],
        order,
    )
    return coarse_on_fine - fine.density[time_index]


def masked_norms(diff: np.ndarray, mask: np.ndarray) -> tuple[float, float]:
    return error_norms(np.where(mask, diff, np.nan))


def build_interp_pairs(runs: list[InterpRun], order: int) -> list[InterpPair]:
    if len(runs) < 2:
        return []

    reference = runs[-1]
    for run in runs[:-1]:
        if len(run.t_norm) != len(reference.t_norm) or not np.allclose(
            run.t_norm,
            reference.t_norm,
            rtol=1.0e-12,
            atol=1.0e-14,
        ):
            raise ValueError(
                f"Time array for {run.case_dir} does not match finest run {reference.case_dir}."
            )

    pairs: list[InterpPair] = []
    for fine_index in range(1, len(runs)):
        coarse = runs[fine_index - 1]
        fine = runs[fine_index]
        compare_mask = interp_overlap_mask(coarse, fine)
        point_count = int(np.count_nonzero(compare_mask))
        if point_count == 0:
            raise ValueError(
                f"No overlapping plasma cells available for {coarse.case_dir} -> {fine.case_dir}."
            )

        l2_series = np.empty(len(fine.t_norm), dtype=float)
        linf_series = np.empty(len(fine.t_norm), dtype=float)
        selected_diff: np.ndarray | None = None

        print(
            f"Interpolating {coarse.case_dir.name} -> {fine.case_dir.name} "
            f"over {point_count} target plasma cells."
        )
        placeholder = InterpPair(
            bc_name=fine.bc_name,
            coarse=coarse,
            fine=fine,
            h=fine.h,
            mask=compare_mask,
            point_count=point_count,
            l2=math.nan,
            linf=math.nan,
            l2_series=l2_series,
            linf_series=linf_series,
            error_selected=np.empty(fine.R.shape, dtype=float),
        )
        for tind in range(len(fine.t_norm)):
            diff = interp_pair_diff(placeholder, tind, order)
            l2_series[tind], linf_series[tind] = masked_norms(diff, compare_mask)
            if tind == fine.time_index:
                selected_diff = diff

        if selected_diff is None:
            selected_diff = interp_pair_diff(placeholder, fine.time_index, order)
        error_selected = np.where(compare_mask, np.abs(selected_diff), np.nan)
        pairs.append(
            InterpPair(
                bc_name=fine.bc_name,
                coarse=coarse,
                fine=fine,
                h=fine.h,
                mask=compare_mask,
                point_count=point_count,
                l2=float(l2_series[fine.time_index]),
                linf=float(linf_series[fine.time_index]),
                l2_series=l2_series,
                linf_series=linf_series,
                error_selected=error_selected,
            )
        )
    return pairs


def fit_interp_order(pairs: list[InterpPair], metric: str = "l2") -> float:
    h = np.asarray([pair.h for pair in pairs], dtype=float)
    err = np.asarray([getattr(pair, metric) for pair in pairs], dtype=float)
    return fit_order_values(h, err)


def selected_interp_pairs(
    neumann: list[InterpPair],
    dirichlet: list[InterpPair],
    use_dirichlet: bool,
) -> list[InterpPair]:
    selected = dirichlet if use_dirichlet else neumann
    if selected:
        return selected
    fallback = neumann if use_dirichlet else dirichlet
    if fallback:
        warnings.warn("Requested interpolation family is missing; plotting the available family.")
        return fallback
    raise FileNotFoundError("No completed interpolation pairs were found.")


def print_interp_summary(neumann: list[InterpPair], dirichlet: list[InterpPair]) -> None:
    print("\nSuccessive self-convergence using target-grid plasma cells")
    print("Interpolation: tensor-product spline on each coarse grid")
    for pairs, label in ((neumann, "Neumann"), (dirichlet, "Dirichlet")):
        if not pairs:
            print(f"\n{label}: no successive pairs")
            continue
        print(f"\n{label}:")
        print("  pair                 N_coarse  N_fine       h_fine        points       L2        Linf")
        for pair in pairs:
            print(
                f"  {pair.label:<18s} {pair.coarse.resolution:8d} {pair.fine.resolution:7d} "
                f"{pair.h:14.8e} {pair.point_count:9d} "
                f"{pair.l2:.6e} {pair.linf:.6e}"
            )
        print(f"  observed L2 order: {fit_interp_order(pairs, 'l2'):.3f}")
        print(f"  observed Linf order: {fit_interp_order(pairs, 'linf'):.3f}")


def plot_interp_summary(
    neumann: list[InterpPair],
    dirichlet: list[InterpPair],
    case_name: str,
    output: Path,
    dpi: int,
    show: bool,
    use_dirichlet: bool,
    log_error: bool,
) -> None:
    plot_pair = selected_interp_pairs(neumann, dirichlet, use_dirichlet)[-1]
    plot_run = plot_pair.fine
    plot_label = bc_display_name(plot_pair.bc_name)
    density = np.where(plot_run.mask, plot_run.density[plot_run.time_index], np.nan)
    error, error_label = error_display_values(plot_pair.error_selected, log_error)

    fig, axes = plt.subplots(1, 3, figsize=(17, 5), constrained_layout=True)

    dens_plot = axes[0].pcolormesh(plot_run.R, plot_run.Z, masked_frame(density), shading="nearest")
    axes[0].set_title(f"{plot_label} plasma density, N={plot_run.resolution}")
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("z")
    axes[0].set_aspect("equal", adjustable="box")
    fig.colorbar(dens_plot, ax=axes[0])

    err_plot = axes[1].pcolormesh(
        plot_run.R,
        plot_run.Z,
        masked_frame(error),
        shading="nearest",
    )
    axes[1].set_title(f"{plot_label} successive {error_label}")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("z")
    axes[1].set_aspect("equal", adjustable="box")
    fig.colorbar(err_plot, ax=axes[1])

    for pairs, label, color, marker in (
        (neumann, "Neumann", "blue", "s"),
        (dirichlet, "Dirichlet", "red", "x"),
    ):
        if not pairs:
            continue
        suffix = label.lower()[:3]
        h = np.asarray([pair.h for pair in pairs], dtype=float)
        l2_err = np.asarray([pair.l2 for pair in pairs], dtype=float)
        linf_err = np.asarray([pair.linf for pair in pairs], dtype=float)
        l2_order = fit_interp_order(pairs, "l2")
        linf_order = fit_interp_order(pairs, "linf")
        axes[2].loglog(
            h,
            l2_err,
            marker=marker,
            color=color,
            linewidth=2,
            label=rf"$L_{{2,{suffix}}}$, p={l2_order:.2f}",
        )
        axes[2].loglog(
            h,
            linf_err,
            marker=marker,
            color=color,
            linestyle="--",
            linewidth=2,
            label=rf"$L_{{\infty,{suffix}}}$, p={linf_order:.2f}",
        )

    axes[2].set_title("Successive convergence, plasma cells")
    axes[2].set_xlabel("Target grid spacing 1/N")
    axes[2].set_ylabel("Successive error")
    axes[2].grid(True, which="both", linestyle="--", alpha=0.4)
    axes[2].legend()

    fig.suptitle(
        f"{case_name} interpolation comparison at t={plot_run.time_seconds:g} "
        f"(normalized {plot_run.time_norm_selected:.3e})"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi)
    if show:
        plt.show()
    else:
        plt.close(fig)


def show_interp_density_animation(
    pair: InterpPair,
    case_name: str,
    order: int,
    log_error: bool,
    save_path: Path | None = None,
    dpi: int = 200,
    show: bool = True,
) -> None:
    fine = pair.fine
    plot_label = bc_display_name(pair.bc_name)
    density = np.where(fine.mask[None, :, :], fine.density, np.nan)
    density_values = density[np.isfinite(density)]
    if not density_values.size:
        raise ValueError(f"No finite plasma density values found for {fine.case_dir}")

    initial_error, error_label = error_display_values(
        np.where(pair.mask, np.abs(interp_pair_diff(pair, 0, order)), np.nan),
        log_error,
    )

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.8))
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.19, top=0.88, wspace=0.28)
    title = fig.suptitle("")

    density_mesh = axes[0].pcolormesh(
        fine.R,
        fine.Z,
        masked_frame(density[0]),
        shading="nearest",
        vmin=float(np.nanmin(density_values)),
        vmax=float(np.nanmax(density_values)),
    )
    axes[0].set_title(f"{plot_label} plasma density, N={fine.resolution}")
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("z")
    axes[0].set_aspect("equal", adjustable="box")
    fig.colorbar(density_mesh, ax=axes[0])

    error_mesh = axes[1].pcolormesh(
        fine.R,
        fine.Z,
        masked_frame(initial_error),
        shading="nearest",
    )
    apply_clim(error_mesh, initial_error)
    axes[1].set_title(f"{plot_label} successive {error_label}")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("z")
    axes[1].set_aspect("equal", adjustable="box")
    fig.colorbar(error_mesh, ax=axes[1])

    slider_ax = fig.add_axes([0.18, 0.065, 0.68, 0.035])
    button_ax = fig.add_axes([0.055, 0.052, 0.09, 0.055])
    frame_slider = Slider(
        slider_ax,
        "Frame",
        0,
        len(fine.t_array) - 1,
        valinit=0,
        valstep=1,
        valfmt="%0.0f",
    )
    play_button = Button(button_ax, "Pause")

    current_frame = 0
    is_playing = True
    updating_slider = False
    animation: FuncAnimation | None = None

    def draw_frame(tind: int, update_slider: bool = True):
        nonlocal current_frame, updating_slider
        tind = int(np.clip(tind, 0, len(fine.t_array) - 1))
        current_frame = tind
        error, _ = error_display_values(
            np.where(pair.mask, np.abs(interp_pair_diff(pair, tind, order)), np.nan),
            log_error,
        )
        density_mesh.set_array(masked_frame(density[tind]).ravel())
        error_mesh.set_array(masked_frame(error).ravel())
        apply_clim(error_mesh, error)
        title.set_text(
            f"{case_name} interpolation animation, {pair.label}, "
            f"frame {tind + 1}/{len(fine.t_array)}, t={fine.t_array[tind]:g} "
            f"(normalized {fine.t_norm[tind]:.3e})"
        )
        if update_slider:
            updating_slider = True
            frame_slider.set_val(tind)
            updating_slider = False
        return density_mesh, error_mesh, title

    def set_playing(value: bool) -> None:
        nonlocal is_playing
        is_playing = value
        play_button.label.set_text("Pause" if is_playing else "Play")
        if animation is not None:
            if is_playing:
                animation.event_source.start()
            else:
                animation.event_source.stop()
        fig.canvas.draw_idle()

    def on_slider(value: float) -> None:
        if updating_slider:
            return
        set_playing(False)
        draw_frame(int(round(value)), update_slider=False)
        fig.canvas.draw_idle()

    def on_play(_event) -> None:
        set_playing(not is_playing)

    def update(_frame):
        return draw_frame((current_frame + 1) % len(fine.t_array))

    draw_frame(0, update_slider=False)
    if save_path is not None:
        save_animation = FuncAnimation(
            fig,
            lambda frame: draw_frame(frame, update_slider=False),
            frames=range(len(fine.t_array)),
            interval=250,
            repeat=False,
            cache_frame_data=False,
        )
        save_animation_file(save_animation, save_path, dpi)
        del save_animation
        draw_frame(0, update_slider=False)

    if not show:
        plt.close(fig)
        return

    frame_slider.on_changed(on_slider)
    play_button.on_clicked(on_play)
    animation = FuncAnimation(
        fig,
        update,
        interval=250,
        repeat=True,
        cache_frame_data=False,
    )
    plt.show()
    plt.close(fig)
    del animation


def nearest_mask_sample(
    mask: np.ndarray,
    R1d: np.ndarray,
    Z1d: np.ndarray,
    RR: np.ndarray,
    ZZ: np.ndarray,
) -> np.ndarray:
    i = np.searchsorted(R1d, RR, side="left")
    k = np.searchsorted(Z1d, ZZ, side="left")
    i = np.clip(i, 1, len(R1d) - 1)
    k = np.clip(k, 1, len(Z1d) - 1)
    i -= np.abs(RR - R1d[i - 1]) <= np.abs(RR - R1d[i])
    k -= np.abs(ZZ - Z1d[k - 1]) <= np.abs(ZZ - Z1d[k])
    return mask[i, k]


def cell_center_plot_fields(run: RunResult) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x_slice = core_x_slice(run)
    mask = run.mask[x_slice, :]
    density = np.where(mask, run.density[x_slice, :], np.nan)
    error = np.where(mask, run.error[x_slice, :], np.nan)
    return run.R[x_slice, :], run.Z[x_slice, :], density, error


def sample_on_regular_grid(run: RunResult, nplot: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x_slice = core_x_slice(run)
    R1d = run.R[x_slice, 0]
    Z1d = run.Z[0, :]
    r = np.linspace(float(np.nanmin(R1d)), float(np.nanmax(R1d)), nplot)
    z = np.linspace(float(np.nanmin(Z1d)), float(np.nanmax(Z1d)), nplot)
    RR, ZZ = np.meshgrid(r, z, indexing="ij")
    density = run.spline_density.ev(RR, ZZ)
    truth = run.spline_truth.ev(RR, ZZ)
    visible = nearest_mask_sample(
        run.mask[x_slice, :],
        R1d,
        Z1d,
        RR,
        ZZ,
    )
    density = np.where(visible, density, np.nan)
    truth = np.where(visible, truth, np.nan)
    return RR, ZZ, density, np.abs(density - truth)


def plot_summary(
    neumann: list[RunResult],
    dirichlet: list[RunResult],
    case_name: str,
    output: Path,
    nplot: int,
    dpi: int,
    show: bool,
    use_dirichlet: bool,
    log_error: bool,
) -> None:
    plot_run = selected_plot_family(neumann, dirichlet, use_dirichlet)[-1]
    plot_label = bc_display_name(plot_run.bc_name)
    R, Z, density, error = cell_center_plot_fields(plot_run)
    error, error_label = error_display_values(error, log_error)

    fig, axes = plt.subplots(1, 3, figsize=(17, 5), constrained_layout=True)

    dens_plot = axes[0].pcolormesh(R, Z, masked_frame(density), shading="nearest")
    axes[0].set_title(f"{plot_label} plasma density, N={plot_run.resolution}")
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("z")
    axes[0].set_aspect("equal", adjustable="box")
    fig.colorbar(dens_plot, ax=axes[0])

    err_plot = axes[1].pcolormesh(R, Z, masked_frame(error), shading="nearest")
    axes[1].set_title(f"{plot_label} plasma {error_label}")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("z")
    axes[1].set_aspect("equal", adjustable="box")
    fig.colorbar(err_plot, ax=axes[1])

    for family, label, color, marker in (
        (neumann, "Neumann", "blue", "s"),
        (dirichlet, "Dirichlet", "red", "x"),
    ):
        suffix = label.lower()[:3]
        h = np.asarray([run.h for run in family], dtype=float)
        l2_err = np.asarray([run.l2 for run in family], dtype=float)
        linf_err = np.asarray([run.linf for run in family], dtype=float)
        l2_order = fit_order(family, "l2")
        linf_order = fit_order(family, "linf")
        axes[2].loglog(
            h,
            l2_err,
            marker=marker,
            color=color,
            linewidth=2,
            label=rf"$L_{{2,{suffix}}}$, p={l2_order:.2f}",
        )
        axes[2].loglog(
            h,
            linf_err,
            marker=marker,
            color=color,
            linestyle="--",
            linewidth=2,
            label=rf"$L_{{\infty,{suffix}}}$, p={linf_order:.2f}",
        )

    axes[2].set_title("Density convergence, plasma cells")
    axes[2].set_xlabel("Grid spacing 1/N")
    axes[2].set_ylabel("Error")
    axes[2].grid(True, which="both", linestyle="--", alpha=0.4)
    axes[2].legend()

    fig.suptitle(
        f"{case_name} MMS at t={plot_run.time_seconds:g} "
        f"(normalized {plot_run.time_norm:.3e})"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi)
    if show:
        plt.show()
    else:
        plt.close(fig)


def masked_frame(values: np.ndarray, mask_nonpositive: bool = False) -> np.ma.MaskedArray:
    values = np.asarray(values, dtype=float)
    mask = ~np.isfinite(values)
    if mask_nonpositive:
        mask |= values <= 0.0
    return np.ma.array(values, mask=mask)


def error_display_values(error: np.ndarray, log_error: bool) -> tuple[np.ndarray, str]:
    error = np.asarray(error, dtype=float)
    if not log_error:
        return error, "density error"

    values = np.full(error.shape, np.nan, dtype=float)
    positive = np.isfinite(error) & (error > 0.0)
    values[positive] = np.log10(error[positive])
    return values, "log10 density error"


def finite_clim(values: np.ndarray) -> tuple[float | None, float | None]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if not finite.size:
        return None, None
    vmin = float(np.nanmin(finite))
    vmax = float(np.nanmax(finite))
    if vmin == vmax:
        pad = 1.0 if vmin == 0.0 else 0.05 * abs(vmin)
        return vmin - pad, vmax + pad
    return vmin, vmax


def apply_clim(mesh, values: np.ndarray) -> None:
    vmin, vmax = finite_clim(values)
    if vmin is not None and vmax is not None:
        mesh.set_clim(vmin, vmax)


def save_animation_file(animation: FuncAnimation, output: Path, dpi: int) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.suffix.lower() == ".gif":
        animation.save(str(output), writer=PillowWriter(fps=4), dpi=dpi)
    else:
        animation.save(str(output), dpi=dpi)
    print(f"Wrote {output}")


def load_animation_fields(run: RunResult):
    bout = parse_bout_inp(run.case_dir / "BOUT.inp")
    D = float(bout["h+"]["anomalous_D"])
    expr = sympy_expression(run.function_expr)
    x_laplacian_scale, z_laplacian_scale = laplacian_coordinate_scales(run.R, run.Z, run.mxg)
    solution_fn, _, _ = build_time_solution(
        expr,
        D,
        run.X,
        run.Zn,
        run.mask,
        x_laplacian_scale,
        z_laplacian_scale,
    )

    density = squeeze_y(collect("Nh+", path=str(run.case_dir), info=False))
    t_array = np.asarray(collect("t_array", path=str(run.case_dir), info=False), dtype=float)
    omega_ci = datafile_scalar(first_dump_file(run.case_dir), "Omega_ci", DEFAULT_OMEGA_CI)
    if omega_ci is None or omega_ci == 0.0:
        omega_ci = DEFAULT_OMEGA_CI
    t_norm = t_array / omega_ci

    x_slice = core_x_slice(run)
    R = run.R[x_slice, :]
    Z = run.Z[x_slice, :]
    X = run.X[x_slice, :]
    Zn = run.Zn[x_slice, :]
    plasma_mask = run.mask[x_slice, :]
    density = np.asarray(density[:, x_slice, :], dtype=float)

    truth = np.empty_like(density)
    for tind, time_value in enumerate(t_norm):
        truth[tind] = grid_values(solution_fn(X, Zn, time_value), R.shape)

    density = np.where(plasma_mask[None, :, :], density, np.nan)
    error = np.where(plasma_mask[None, :, :], np.abs(density - truth), np.nan)
    return R, Z, density, error, t_array, t_norm


def show_density_animation(
    run: RunResult,
    case_name: str,
    log_error: bool,
    save_path: Path | None = None,
    dpi: int = 200,
    show: bool = True,
) -> None:
    R, Z, density, error, t_array, t_norm = load_animation_fields(run)
    plot_label = bc_display_name(run.bc_name)

    density_values = density[np.isfinite(density)]
    if not density_values.size:
        raise ValueError(f"No finite plasma density values found for {run.case_dir}")

    error_plot, error_label = error_display_values(error, log_error)
    err_vmin, err_vmax = finite_clim(error_plot)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.8))
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.19, top=0.88, wspace=0.28)
    title = fig.suptitle("")

    density_mesh = axes[0].pcolormesh(
        R,
        Z,
        masked_frame(density[0]),
        shading="nearest",
        vmin=float(np.nanmin(density_values)),
        vmax=float(np.nanmax(density_values)),
    )
    axes[0].set_title(f"{plot_label} plasma density, N={run.resolution}")
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("z")
    axes[0].set_aspect("equal", adjustable="box")
    fig.colorbar(density_mesh, ax=axes[0])

    error_mesh = axes[1].pcolormesh(
        R,
        Z,
        masked_frame(error_plot[0]),
        shading="nearest",
        vmin=err_vmin,
        vmax=err_vmax,
    )
    axes[1].set_title(f"{plot_label} plasma {error_label}")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("z")
    axes[1].set_aspect("equal", adjustable="box")
    fig.colorbar(error_mesh, ax=axes[1])

    slider_ax = fig.add_axes([0.18, 0.065, 0.68, 0.035])
    button_ax = fig.add_axes([0.055, 0.052, 0.09, 0.055])
    frame_slider = Slider(
        slider_ax,
        "Frame",
        0,
        len(t_array) - 1,
        valinit=0,
        valstep=1,
        valfmt="%0.0f",
    )
    play_button = Button(button_ax, "Pause")

    current_frame = 0
    is_playing = True
    updating_slider = False
    animation: FuncAnimation | None = None

    def draw_frame(tind: int, update_slider: bool = True):
        nonlocal current_frame, updating_slider
        tind = int(np.clip(tind, 0, len(t_array) - 1))
        current_frame = tind
        density_mesh.set_array(masked_frame(density[tind]).ravel())
        error_mesh.set_array(masked_frame(error_plot[tind]).ravel())
        title.set_text(
            f"{case_name} MMS animation, frame {tind + 1}/{len(t_array)}, "
            f"t={t_array[tind]:g} (normalized {t_norm[tind]:.3e})"
        )
        if update_slider:
            updating_slider = True
            frame_slider.set_val(tind)
            updating_slider = False
        return density_mesh, error_mesh, title

    def set_playing(value: bool) -> None:
        nonlocal is_playing
        is_playing = value
        play_button.label.set_text("Pause" if is_playing else "Play")
        if animation is not None:
            if is_playing:
                animation.event_source.start()
            else:
                animation.event_source.stop()
        fig.canvas.draw_idle()

    def on_slider(value: float) -> None:
        if updating_slider:
            return
        set_playing(False)
        draw_frame(int(round(value)), update_slider=False)
        fig.canvas.draw_idle()

    def on_play(_event) -> None:
        set_playing(not is_playing)

    def update(_frame):
        return draw_frame((current_frame + 1) % len(t_array))

    draw_frame(0, update_slider=False)
    if save_path is not None:
        save_animation = FuncAnimation(
            fig,
            lambda frame: draw_frame(frame, update_slider=False),
            frames=range(len(t_array)),
            interval=250,
            repeat=False,
            cache_frame_data=False,
        )
        save_animation_file(save_animation, save_path, dpi)
        del save_animation
        draw_frame(0, update_slider=False)

    if not show:
        plt.close(fig)
        return

    frame_slider.on_changed(on_slider)
    play_button.on_clicked(on_play)
    animation = FuncAnimation(
        fig,
        update,
        interval=250,
        repeat=True,
        cache_frame_data=False,
    )
    plt.show()
    plt.close(fig)
    del animation


def point_array(grid: dict[str, np.ndarray | float | int], name: str) -> np.ndarray:
    if name not in grid:
        return np.empty((0, 2), dtype=float)
    pts = np.asarray(grid[name], dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 2:
        return np.empty((0, 2), dtype=float)
    return pts[np.all(np.isfinite(pts), axis=1)]


def choose_setup_zoom(
    x_edges: np.ndarray,
    z_edges: np.ndarray,
    bndry_pts: np.ndarray,
    ghost_pts: np.ndarray,
) -> tuple[tuple[float, float], tuple[float, float]]:
    dx = float(np.nanmedian(np.diff(x_edges)))
    dz = float(np.nanmedian(np.diff(z_edges)))
    half_width = 7.5 * max(abs(dx), abs(dz))

    if bndry_pts.size:
        target = np.array([np.nanmedian(bndry_pts[:, 0]), np.nanmax(bndry_pts[:, 1])])
        anchor = bndry_pts[np.nanargmin(np.sum((bndry_pts - target) ** 2, axis=1))]
    elif ghost_pts.size:
        target = np.array([np.nanmedian(ghost_pts[:, 0]), np.nanmax(ghost_pts[:, 1])])
        anchor = ghost_pts[np.nanargmin(np.sum((ghost_pts - target) ** 2, axis=1))]
    else:
        anchor = np.array(
            [
                0.5 * (float(x_edges[0]) + float(x_edges[-1])),
                0.5 * (float(z_edges[0]) + float(z_edges[-1])),
            ]
        )

    xlim = (float(anchor[0] - half_width), float(anchor[0] + half_width))
    ylim = (float(anchor[1] - half_width), float(anchor[1] + half_width))
    xlim = (max(float(x_edges[0]), xlim[0]), min(float(x_edges[-1]), xlim[1]))
    ylim = (max(float(z_edges[0]), ylim[0]), min(float(z_edges[-1]), ylim[1]))
    return xlim, ylim


def points_in_view(points: np.ndarray, xlim: tuple[float, float], ylim: tuple[float, float]) -> bool:
    if points.size == 0:
        return False
    return bool(
        np.any(
            (points[:, 0] >= xlim[0])
            & (points[:, 0] <= xlim[1])
            & (points[:, 1] >= ylim[0])
            & (points[:, 1] <= ylim[1])
        )
    )


def plot_grid_setup(ax, run: RunResult) -> None:
    grid = read_grid(run.grid_path)
    R = squeeze_grid(np.asarray(grid["R"], dtype=float))
    Z = squeeze_grid(np.asarray(grid["Z"], dtype=float))
    x_slice = core_x_slice(run)
    R_core = R[x_slice, :]
    Z_core = Z[x_slice, :]
    x_edges = centers_to_edges(R_core[:, 0])
    z_edges = centers_to_edges(Z_core[0, :])

    in_mask = squeeze_grid(np.asarray(grid["in_mask"])).astype(bool)[x_slice, :]
    ghost_id = squeeze_grid(np.asarray(grid["ghost_id"]))[x_slice, :] if "ghost_id" in grid else None
    ghost_mask = ghost_id >= 0 if ghost_id is not None else np.zeros_like(in_mask, dtype=bool)
    if "cc_mask" in grid:
        cut_mask = squeeze_grid(np.asarray(grid["cc_mask"])).astype(bool)[x_slice, :]
    else:
        vol_frac = squeeze_grid(np.asarray(grid["vol_frac"], dtype=float))[x_slice, :]
        cut_mask = (vol_frac > 0.0) & (vol_frac < 1.0)

    layout = np.full(R_core.shape, np.nan, dtype=float)
    layout[in_mask] = 0.0
    layout[cut_mask] = 1.0
    layout[ghost_mask] = 2.0

    cmap = ListedColormap(["#d8ecff", "#f3b562", "#55a868"])
    cmap.set_bad((1.0, 1.0, 1.0, 0.0))
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)

    ax.vlines(x_edges, z_edges[0], z_edges[-1], color="0.86", linewidth=0.45, zorder=0)
    ax.hlines(z_edges, x_edges[0], x_edges[-1], color="0.86", linewidth=0.45, zorder=0)
    ax.pcolormesh(
        R_core,
        Z_core,
        np.ma.masked_invalid(layout),
        shading="nearest",
        cmap=cmap,
        norm=norm,
        zorder=1,
    )

    ghost_pts = point_array(grid, "ghost_pts")
    image_pts = point_array(grid, "image_pts")
    bndry_pts = point_array(grid, "bndry_pts")
    xlim, ylim = choose_setup_zoom(x_edges, z_edges, bndry_pts, ghost_pts)

    boundary = point_array(grid, "eff_bndry_vertices")
    if boundary.size:
        ax.plot(boundary[:, 0], boundary[:, 1], color="black", linewidth=1.6, zorder=4)
    elif bndry_pts.size:
        ax.plot(bndry_pts[:, 0], bndry_pts[:, 1], color="black", linewidth=1.2, zorder=4)

    if ghost_pts.size and image_pts.shape == ghost_pts.shape and bndry_pts.shape == ghost_pts.shape:
        for ghost, boundary_pt, image in zip(ghost_pts, bndry_pts, image_pts):
            pts = np.vstack((ghost, boundary_pt, image))
            if points_in_view(pts, xlim, ylim):
                ax.plot(pts[:, 0], pts[:, 1], color="crimson", linewidth=0.8, alpha=0.45, zorder=3)

    if ghost_pts.size:
        ax.scatter(ghost_pts[:, 0], ghost_pts[:, 1], s=14, c="#267f3e", marker="o", zorder=5)
    if image_pts.size:
        ax.scatter(image_pts[:, 0], image_pts[:, 1], s=12, c="crimson", marker="x", zorder=5)
    if bndry_pts.size:
        ax.scatter(bndry_pts[:, 0], bndry_pts[:, 1], s=10, c="black", marker=".", zorder=5)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(f"IB setup, N={run.resolution}")
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    ax.legend(
        handles=[
            Patch(facecolor="#d8ecff", edgecolor="none", label="plasma"),
            Patch(facecolor="#f3b562", edgecolor="none", label="cut cell"),
            Patch(facecolor="#55a868", edgecolor="none", label="IB ghost"),
            Line2D([0], [0], color="black", linewidth=1.6, label="boundary"),
            Line2D([0], [0], color="crimson", linewidth=0.8, label="ghost map"),
        ],
        fontsize=7,
        loc="best",
    )


def plot_ghost_summary(
    neumann: list[RunResult],
    dirichlet: list[RunResult],
    case_name: str,
    output: Path,
    dpi: int,
    show: bool,
    use_dirichlet: bool,
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True)
    panels = (
        ("t=0", lambda run: run.ghost_t0),
        ("selected time", lambda run: run.ghost_selected),
    )

    for ax, (title, getter) in zip(axes, panels):
        for family, label, color, marker in (
            (neumann, "Neumann", "blue", "s"),
            (dirichlet, "Dirichlet", "red", "x"),
        ):
            suffix = label.lower()[:3]
            h = np.asarray([run.h for run in family], dtype=float)
            l2_err = np.asarray([getter(run).l2 for run in family], dtype=float)
            linf_err = np.asarray([getter(run).linf for run in family], dtype=float)
            l2_order = fit_order_values(h, l2_err)
            linf_order = fit_order_values(h, linf_err)
            ax.loglog(
                h,
                l2_err,
                marker=marker,
                color=color,
                linewidth=2,
                label=rf"$L_{{2,{suffix}}}$, p={l2_order:.2f}",
            )
            ax.loglog(
                h,
                linf_err,
                marker=marker,
                color=color,
                linestyle="--",
                linewidth=2,
                label=rf"$L_{{\infty,{suffix}}}$, p={linf_order:.2f}",
            )

        ax.set_title(f"Ghost value convergence, {title}")
        ax.set_xlabel("Grid spacing 1/N")
        ax.set_ylabel("Error")
        ax.grid(True, which="both", linestyle="--", alpha=0.4)
        ax.legend()

    setup_run = selected_plot_family(neumann, dirichlet, use_dirichlet)[0]
    plot_grid_setup(axes[2], setup_run)

    selected = selected_plot_family(neumann, dirichlet, use_dirichlet)[-1].ghost_selected
    fig.suptitle(
        f"{case_name} analytic ghost-value MMS; selected t={selected.time_seconds:g} "
        f"(normalized {selected.time_norm:.3e})"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi)
    if show:
        plt.show()
    else:
        plt.close(fig)


def print_summary(neumann: list[RunResult], dirichlet: list[RunResult]) -> None:
    for family, label in ((neumann, "Neumann"), (dirichlet, "Dirichlet")):
        print(f"\n{label}:")
        for run in family:
            lam = "none" if run.lambda_value is None else f"{run.lambda_value:.12g}"
            print(
                f"  N={run.resolution:5d} h={run.h:.6e} "
                f"plasma L2={run.l2:.6e} plasma Linf={run.linf:.6e} "
                f"lambda={lam} model={run.time_model}"
            )
        print(f"  observed plasma L2 order: {fit_order(family, 'l2'):.3f}")
        print(f"  observed plasma Linf order: {fit_order(family, 'linf'):.3f}")


def print_ghost_summary(neumann: list[RunResult], dirichlet: list[RunResult]) -> None:
    for getter, label in (
        (lambda run: run.ghost_t0, "Ghosts at t=0"),
        (lambda run: run.ghost_selected, "Ghosts at selected time"),
    ):
        print(f"\n{label}:")
        for family, family_label in ((neumann, "Neumann"), (dirichlet, "Dirichlet")):
            h = np.asarray([run.h for run in family], dtype=float)
            l2_err = np.asarray([getter(run).l2 for run in family], dtype=float)
            linf_err = np.asarray([getter(run).linf for run in family], dtype=float)
            n_ghost = getter(family[-1]).n_ghost
            print(
                f"  {family_label}: ghost L2 order={fit_order_values(h, l2_err):.3f} "
                f"ghost Linf order={fit_order_values(h, linf_err):.3f} "
                f"n_ghost(finest)={n_ghost}"
            )


def output_path(path: Path | None, default_path: Path) -> Path:
    if path is None:
        path = default_path
    if path.is_absolute():
        return path
    return Path.cwd() / path


def animation_output_path(value: str | None, default_path: Path) -> Path | None:
    if value is None:
        return None
    if value == SAVE_ANIM_DEFAULT:
        return output_path(None, default_path)
    return output_path(Path(value), default_path)


def run_interp_mode(
    args: argparse.Namespace,
    base: Path,
    output: Path,
    anim_output: Path | None,
) -> None:
    neumann_runs = load_interp_family(base, "neu", args.time_index)
    dirichlet_runs = load_interp_family(base, "dir", args.time_index)
    neumann = build_interp_pairs(neumann_runs, args.interp_order)
    dirichlet = build_interp_pairs(dirichlet_runs, args.interp_order)
    if not neumann and not dirichlet:
        raise FileNotFoundError(f"No completed interpolation pairs found under {base}")

    print_interp_summary(neumann, dirichlet)
    plot_interp_summary(
        neumann,
        dirichlet,
        base.name,
        output,
        args.dpi,
        args.plot,
        args.plot_dirichlet,
        args.log_error,
    )
    print(f"\nWrote {output}")
    if args.ghost_output is not None:
        warnings.warn("--ghost-output is ignored in --interp mode; no analytic ghost truth is available.")
    if args.anim or anim_output is not None:
        show_interp_density_animation(
            selected_interp_pairs(neumann, dirichlet, args.plot_dirichlet)[-1],
            base.name,
            args.interp_order,
            args.log_error,
            anim_output,
            args.dpi,
            args.anim,
        )


def main() -> None:
    args = parse_args()
    directory = args.directory_option or args.directory or "diamond"
    base = resolve_case_dir(directory)
    mode_suffix = "Interp" if args.interp else ""
    output = output_path(args.output, args.out_dir / f"{base.name}{mode_suffix}MMS.pdf")
    ghost_output = output_path(args.ghost_output, args.out_dir / f"{base.name}GhostMMS.pdf")
    anim_default_name = f"{base.name}{mode_suffix}Anim.gif"
    anim_output = animation_output_path(args.save_anim, args.out_dir / anim_default_name)

    if args.interp:
        run_interp_mode(args, base, output, anim_output)
        return

    neumann = load_family(base, "neu", args.time_index)
    dirichlet = load_family(base, "dir", args.time_index)
    print_summary(neumann, dirichlet)
    print_ghost_summary(neumann, dirichlet)
    plot_summary(
        neumann,
        dirichlet,
        base.name,
        output,
        args.plot_points,
        args.dpi,
        args.plot,
        args.plot_dirichlet,
        args.log_error,
    )
    plot_ghost_summary(
        neumann,
        dirichlet,
        base.name,
        ghost_output,
        args.dpi,
        args.plot,
        args.plot_dirichlet,
    )
    print(f"\nWrote {output}")
    print(f"Wrote {ghost_output}")
    if args.anim or anim_output is not None:
        show_density_animation(
            selected_plot_family(neumann, dirichlet, args.plot_dirichlet)[-1],
            base.name,
            args.log_error,
            anim_output,
            args.dpi,
            args.anim,
        )


if __name__ == "__main__":
    main()
