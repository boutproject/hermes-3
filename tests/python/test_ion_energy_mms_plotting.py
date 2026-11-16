"""Publication artifacts preserve domain and interior convergence measurements."""

import csv
import sys
from pathlib import Path

import pytest

pytest.importorskip("matplotlib")
HERE = Path(__file__).resolve().parents[1] / "integrated/ion-energy-mms"
sys.path.insert(0, str(HERE))
from plotting import plot_convergence


def test_vector_figures_and_csv_preserve_wall_degradation(tmp_path, monkeypatch):
    from matplotlib.figure import Figure

    original_save = Figure.savefig
    ticks = []

    def save(figure, *args, **kwargs):
        original_save(figure, *args, **kwargs)
        for axis in figure.axes:
            ticks.append([label.get_text() for label in axis.get_xticklabels()])
            assert axis.get_xticklabels(minor=True) == []

    monkeypatch.setattr(Figure, "savefig", save)
    errors = {
        size: {
            "Vort": {
                "l2": 0.01 * (8 / size),
                "linf": 0.02,
                "l2_interior": 0.005 * (8 / size) ** 2,
                "linf_interior": 0.01 * (8 / size) ** 2,
            },
            "ion_chi": {
                "l2": 0.02 * (8 / size) ** 2,
                "linf": 0.04 * (8 / size) ** 2,
                "l2_interior": 0.01 * (8 / size) ** 2,
                "linf_interior": 0.02 * (8 / size) ** 2,
                "reference_rms": 0.01,
                "relative_rms_error": 2 * (8 / size) ** 2,
            },
        }
        for size in (8, 16, 32, 64)
    }
    artifacts = plot_convergence(tmp_path, "parallel_current_wall", errors)
    assert all(labels == ["8", "16", "32", "64"] for labels in ticks)
    assert set(artifacts) == {
        "convergence.csv",
        "convergence-Vort.pdf",
        "convergence-Vort.svg",
        "convergence-ion_chi.pdf",
        "convergence-ion_chi.svg",
    }
    with (tmp_path / "convergence.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 8
    assert {row["variant"] for row in rows} == {"parallel_current_wall"}
    finest = next(
        row for row in rows if row["field"] == "Vort" and row["cells"] == "64"
    )
    assert float(finest["linf"]) == 0.02
    assert float(finest["l2_interior"]) == pytest.approx(0.005 / 64)
    assert float(finest["order_l2"]) == 1
    assert float(finest["order_linf"]) == 0
    assert float(finest["order_l2_interior"]) == 2
    for name in ("Vort", "ion_chi"):
        assert (tmp_path / f"convergence-{name}.pdf").read_bytes().startswith(b"%PDF")
        svg = (tmp_path / f"convergence-{name}.svg").read_text()
        assert "<svg" in svg and "<text" in svg
        assert "Domain" in svg and "Interior" in svg


def test_partial_and_zero_error_output(tmp_path):
    errors = {
        8: {"phi": dict.fromkeys(("l2", "linf", "l2_interior", "linf_interior"), 0)}
    }
    artifacts = plot_convergence(tmp_path, "viscosity", errors)
    assert len(artifacts) == 3
    with (tmp_path / "convergence.csv").open() as stream:
        [row] = csv.DictReader(stream)
    assert row["order_l2"] == ""


@pytest.mark.parametrize("sizes, ratio", [((8, 16), 0.4), ((8,), 0.999)])
def test_narrow_error_ranges_retain_numeric_y_labels(
    tmp_path, monkeypatch, sizes, ratio
):
    from matplotlib.figure import Figure

    original_save = Figure.savefig

    def save(figure, *args, **kwargs):
        original_save(figure, *args, **kwargs)
        for axis in figure.axes:
            lower, upper = axis.get_ylim()
            labeled = [
                tick
                for tick, label in zip(axis.get_yticks(), axis.get_yticklabels())
                if lower <= tick <= upper and label.get_text()
            ]
            assert len(labeled) >= 2
            assert all(tick > 0 for tick in labeled)

    monkeypatch.setattr(Figure, "savefig", save)
    errors = {
        size: {
            "Vort": {
                norm: 0.0137 * (ratio if size == 16 or "interior" in norm else 1)
                for norm in ("l2", "linf", "l2_interior", "linf_interior")
            }
        }
        for size in sizes
    }
    assert len(plot_convergence(tmp_path, "parallel_current", errors)) == 3
