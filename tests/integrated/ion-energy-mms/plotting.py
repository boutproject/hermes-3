"""Publication convergence figures from recorded MMS norms; no simulation reads."""

import csv

import numpy as np


def plot_convergence(output, variant, errors):
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.ticker import (
        FixedLocator,
        LogFormatterSciNotation,
        LogLocator,
        MaxNLocator,
        NullFormatter,
        NullLocator,
    )

    sizes = sorted(errors)
    fields = list(errors[sizes[0]])
    artifacts = []
    norms = ("l2", "linf", "l2_interior", "linf_interior")
    table = output / "convergence.csv"
    with table.open("w", newline="") as stream:
        columns = [
            "variant",
            "cells",
            "field",
            "l2",
            "linf",
            "l2_interior",
            "linf_interior",
            "reference_rms",
            "relative_rms_error",
            *("order_" + norm for norm in norms),
        ]
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for index, size in enumerate(sizes):
            for name in fields:
                orders = {}
                for norm in norms:
                    order = None
                    if index:
                        coarse = errors[sizes[index - 1]][name][norm]
                        fine = errors[size][name][norm]
                        if (
                            np.isfinite([coarse, fine]).all()
                            and coarse > 0
                            and fine > 0
                        ):
                            order = float(
                                np.log(coarse / fine) / np.log(size / sizes[index - 1])
                            )
                    orders["order_" + norm] = order
                writer.writerow(
                    {
                        "variant": variant,
                        "cells": size,
                        "field": name,
                        **{k: errors[size][name].get(k) for k in columns[3:]},
                        **orders,
                    }
                )
    artifacts.append(table.name)
    with plt.rc_context({"font.size": 10, "svg.fonttype": "none", "pdf.fonttype": 42}):
        for name in fields:
            fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.8), layout="constrained")
            for axis, norm, label in zip(
                axes, ("l2", "linf"), (r"$L_2$", r"$L_\infty$")
            ):
                plotted = False
                for suffix, region, marker in (
                    ("", "Domain", "o"),
                    ("_interior", "Interior", "s"),
                ):
                    values = np.array(
                        [errors[size][name][norm + suffix] for size in sizes]
                    )
                    valid = np.isfinite(values) & (values > 0)
                    plotted |= bool(valid.any())
                    axis.plot(
                        np.array(sizes)[valid],
                        values[valid],
                        marker=marker,
                        label=region,
                    )
                first = errors[sizes[0]][name][norm]
                if first > 0 and len(sizes) > 1:
                    for order, style in ((1, ":"), (2, "--")):
                        axis.plot(
                            sizes,
                            first * (np.array(sizes) / sizes[0]) ** (-order),
                            style,
                            color="0.6",
                            linewidth=1,
                            label=f"Order {order}",
                        )
                axis.set(
                    xscale="log",
                    yscale="log",
                    xlabel="X/Z cells",
                    ylabel=label + " error",
                )
                if not plotted:
                    # Exact-zero legacy fields and partial reports have no
                    # positive data from which matplotlib can infer log limits.
                    axis.set_ylim(0.1, 10)
                    axis.set_yticks([])
                    axis.text(
                        0.5,
                        0.5,
                        "No positive error",
                        ha="center",
                        transform=axis.transAxes,
                    )
                axis.set_xticks(sizes, labels=[str(size) for size in sizes])
                axis.xaxis.set_minor_locator(NullLocator())
                axis.yaxis.set_minor_formatter(NullFormatter())
                if plotted:
                    lower, upper = axis.get_ylim()
                    ticks = axis.yaxis.get_major_locator().tick_values(lower, upper)
                    if np.count_nonzero((ticks >= lower) & (ticks <= upper)) < 2:
                        # Two-level studies can span less than a decade. Label
                        # intermediate values so these panels keep a numeric scale.
                        ticks = LogLocator(subs=(1, 2, 5)).tick_values(lower, upper)
                        ticks = ticks[(ticks >= lower) & (ticks <= upper)]
                        if len(ticks) < 2:
                            ticks = MaxNLocator(nbins=3).tick_values(lower, upper)
                            ticks = ticks[(ticks >= lower) & (ticks <= upper)]
                        axis.yaxis.set_major_locator(FixedLocator(ticks))
                        axis.yaxis.set_major_formatter(
                            LogFormatterSciNotation(minor_thresholds=(np.inf, np.inf))
                        )
                axis.grid(visible=True, which="major", alpha=0.2)
            axes[0].legend(frameon=False, fontsize=8)
            fig.suptitle(name + ": " + variant.replace("_", " "))
            for extension in ("pdf", "svg"):
                path = output / f"convergence-{name}.{extension}"
                fig.savefig(path, metadata={"Title": name + " convergence: " + variant})
                artifacts.append(path.name)
            plt.close(fig)
    return artifacts
