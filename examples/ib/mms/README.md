# IB MMS workflow

This folder contains a small workflow for generating immersed-boundary MMS grids,
running Hermes-3 cases, and plotting convergence diagnostics. The focused tools
are `grid_gen.py`, `run.py`, and `mms.py`; `run_all.py` ties them together for a
full sweep.

## 1. Generate grids

Use `grid_gen.py` to generate Zoidberg FCI grid files under
`<case>/gridfiles`.

```bash
python -B grid_gen.py --d circle --np 64 --n 4
python -B grid_gen.py --d diamond --np 64 --n 4
```

Important flags:

- `--d`, `--case`: MMS geometry/case to generate. Current choices are `circle` and
  `diamond`.
- `--np`: starting perpendicular resolution in `x` and `z`. For example,
  `--np 64`.
- `--n`: number of resolution levels, doubling from `--np`. For example,
  `--np 64 --n 4` generates `64`, `128`, `256`, and `512`.
- `--force`: overwrite existing grid files.
- `--dry-run`: print the Zoidberg calls without generating files.

For `circle`, both the centered and shifted circle grids are generated:
`cl<N>.fci.nc` and `clShift<N>.fci.nc`. For `diamond`, the generated files are
`sqrt<N>.fci.nc`.

## 2. Run cases

Use `run.py` to launch Hermes-3 on case directories whose `BOUT.inp` references
an existing grid file.

```bash
python -B run.py --d circle -n 8
python -B run.py --d diamond -n 8
python -B run.py --d circle --dir -n 8
```

The command run for each case is:

```bash
mpirun -n <nproc> --use-hwthread-cpus <git_root>/build/hermes-3 -d <run_dir>
```

On Perlmutter, `run.py` uses `srun` by default instead:

```bash
srun -n <nproc> <git_root>/build/hermes-3 -d <run_dir>
```

Important flags:

- `directory`: optional positional case directory or case name under this
  folder. Defaults to `diamond`.
- `--d`, `--case`, `--case-dir`: case directory or case name under this folder.
  Overrides the positional `directory`.
- `-n`, `--nproc`: number of MPI processes. If omitted, `run.py` chooses a safe
  value from each gridfile, capped by the detected Slurm allocation or local CPU
  count.
- `--launcher`: choose `auto`, `mpirun`, or `srun`. `auto` uses `srun` on
  Perlmutter and `mpirun` otherwise.
- `--auto-nproc`: choose or cap the number of MPI processes from each gridfile.
  The cap is based on `(nx - 2*MXG) / (2*MXG)`, adjusted downward so the number
  of ranks divides the physical x cells. Automatic sizing is also capped by the
  detected Slurm allocation or local CPU count.
- `--max-nproc`: override the detected maximum MPI process count used by
  automatic sizing. For example, use `--max-nproc 128` to cap every generated
  `srun -n` command at 128 ranks.
- `--dir`: run only Dirichlet cases.
- `--neu`: run only Neumann cases.
- `--no-prepare`: do not generate/update resolution run directories from the
  `dir/BOUT.inp` and `neu/BOUT.inp` templates.
- `--prepare-only`: generate/update resolution run directories, then exit.
- `--dry-run`: print the commands without running Hermes.

`run.py` first uses the `dir/BOUT.inp` and `neu/BOUT.inp` files as templates. It
infers matching grid files from each template's `[mesh] file` entry, creates
resolution folders such as `dir/64` or `neu/512`, and writes the per-resolution
`BOUT.inp` with the correct relative grid path. It then searches for leaf
directories with `BOUT.inp`, skips cases with missing grids, and stops if a
launched case returns a nonzero exit code.

## 3. Plot MMS diagnostics

Use `mms.py` after the runs have produced `BOUT.dmp.*.nc` files.

```bash
python -B mms.py --d circle
python -B mms.py --d diamond --plot
python -B mms.py --d diamond --anim
python -B mms.py --d diamond --save-anim
```

Default outputs go into `./out`:

```text
out/<case>MMS.pdf
out/<case>GhostMMS.pdf
out/<case>InterpMMS.pdf  # with --interp
```

Important flags:

- `--d`, `--case`, `--case-dir`: case directory or case name, such as `circle`
  or `diamond`.
- `--tind`: time index to analyze. Defaults to the halfway timestep, `nt // 2`.
  Negative indices work, so `--tind -1` selects the final stored time.
- `--out-dir`: directory for default PDFs. Defaults to `./out`.
- `--output`: explicit density/convergence PDF path.
- `--ghost-output`: explicit ghost convergence/setup PDF path.
- `--plot`, `--show`: show the static figures interactively after saving.
- `--anim`: after saving the static PDFs, show an interactive animation of the
  highest-resolution density and density error over all stored times.
- `--save-anim [PATH]`: save the highest-resolution density/error animation.
  With no path, this writes `out/<case>Anim.gif`, or
  `out/<case>InterpAnim.gif` in `--interp` mode. Combine with `--anim` to save
  and then show the interactive window.
- `--dir`: use Dirichlet cases for the single-family density, error, setup, and
  animation panels. The convergence curves still show both Neumann and
  Dirichlet.
- `--interp`: compare each run against the next higher-resolution run instead
  of using an analytic MMS solution. This is intended for non-analytic cases
  such as `d3d`; it plots successive self-convergence over target-grid plasma
  cells and does not generate an analytic ghost-value PDF.
- `--interp-order`: spline order for `--interp` comparisons. Defaults to `1`,
  matching the old notebook's bilinear comparison.
- `--log`: plot spatial density error panels and error animations as
  `log10(error)`. Nonpositive error values are hidden in these panels.
- `--dpi`: PDF raster DPI. Defaults to `200`.
- `--plot-points`: deprecated compatibility option. The analytic density/error
  panels now plot raw masked grid cells so ghost/outside values cannot pollute
  the boundary through spline resampling.

The main MMS density/error calculations and spatial density/error panels use
plasma cells only: BOUT x-guards and IB ghost cells are excluded. The ghost
convergence plot is separate and compares stored ghost values against the
analytic MMS value at the ghost points.

The reported/plotted norms are:

- `L2`: RMS error, `sqrt(mean(error**2))`, over the selected cells.
- `Linf`: maximum absolute error, `max(abs(error))`.

## Common examples

Generate circle grids from `64` through `512`, run them on 8 MPI ranks, and make
plots:

```bash
python -B grid_gen.py --d circle --np 64 --n 4
python -B run.py --d circle -n 8
python -B mms.py --d circle
```

Preview a run sweep without launching Hermes:

```bash
python -B run.py --d diamond -n 8 --dry-run
```

Preview the Perlmutter launcher and grid-based rank selection:

```bash
python -B run.py --d diamond --launcher srun --auto-nproc --dry-run
```

Cap automatic Perlmutter rank selection to a known allocation size:

```bash
python -B run.py --d diamond --launcher srun --auto-nproc --max-nproc 128 --dry-run
```

Rerun only the circle Dirichlet cases:

```bash
python -B run.py --d circle --dir -n 8
```

Plot the final timestep and show the figures:

```bash
python -B mms.py --d diamond --tind -1 --plot
```

Use Dirichlet data for the single-family spatial plots and animation:

```bash
python -B mms.py --d diamond --dir --anim
```

Run interpolation/self-convergence analysis for a non-analytic case:

```bash
python -B mms.py --d d3d --interp
python -B mms.py --d d3d --interp --dir --tind -1
python -B mms.py --d diamond --interp --log
```

## 4. Run the full workflow

Use `run_all.py` to generate grids, prepare/run Hermes cases, and generate both
Neumann-selected and Dirichlet-selected MMS plots.

```bash
python -B run_all.py -n 8
python -B run_all.py --d circle -n 8 --np 64 --levels 4
```

By default, `run_all.py` processes every case known to `grid_gen.py` plus case
directories that contain `dir` or `neu` subdirectories. Local default sweeps
skip the `d3d` skeleton; on Perlmutter it is included. Use `--d d3d` to request
it explicitly on any system.

Important flags:

- `--d`, `--case`: process one case. Repeat the flag for multiple cases.
- `--np`: starting perpendicular gridpoints for `grid_gen.py`. Defaults to `64`,
  except `d3d` defaults to `256`.
- `--levels`, `--nlevels`, `--grid-n`: number of resolution levels for
  `grid_gen.py`. Defaults to `4`, or `6` on Perlmutter for regular cases such
  as `circle` and `diamond`. `d3d` defaults to `4`.
- `-n`, `--nproc`: number of MPI processes for `run.py`. Defaults to `run.py`'s
  host-aware choice.
- `--launcher`: launcher passed to `run.py`: `auto`, `mpirun`, or `srun`.
- `--auto-nproc`: pass `--auto-nproc` to `run.py`.
- `--max-nproc`: maximum MPI processes passed to `run.py` automatic sizing.
- `--tind`: time index passed to `mms.py`.
- `--out-dir`: directory for generated PDFs. Defaults to `./out`.
- `--interp`: pass `--interp` to `mms.py` for every selected case. `d3d` uses
  interpolation mode by default.
- `--log`: pass `--log` to `mms.py`.
- `--force-grids`: regenerate existing grids.
- `--skip-grid`, `--skip-run`, `--skip-plot`: skip individual workflow stages.
- `--plot`, `--show`: show `mms.py` plots interactively after saving.
- `--dry-run`: print the commands without running them.

`d3d` grid generation is skipped until it is configured in `grid_gen.py`, but
`run_all.py --d d3d --dry-run` will still show the run and interpolation-plot
commands. Once the d3d grid generator is registered, `run_all.py --d d3d` will
use `--np 256 --levels 4` by default. In non-dry-run mode, unconfigured case
directories with no `BOUT.inp` files or no output dumps are skipped cleanly.

The plot outputs are named so the Neumann-selected and Dirichlet-selected
diagnostics do not overwrite each other:

```text
out/<case>NeuMMS.pdf
out/<case>NeuGhostMMS.pdf
out/<case>DirMMS.pdf
out/<case>DirGhostMMS.pdf
```

For interpolation-mode cases such as `d3d`, the run-all outputs include
`Interp` in the name, such as `out/<case>NeuInterpMMS.pdf` and
`out/<case>DirInterpMMS.pdf`. No ghost PDF is written because there is no
analytic ghost-value reference.
