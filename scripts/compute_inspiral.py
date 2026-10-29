#!/usr/bin/env python3
"""Evolve the two-component inspiral model and write trace/summary results.

Run with ``.venv/bin/python scripts/compute_inspiral.py``; default paths resolve
from the repository root. The script always reads the saved APR background
(``results/summary.json`` and ``results/profiles.csv``), so run
``scripts/compute_structure.py`` first. Outputs in ``--output`` (default
``results/inspiral/``): one CSV trace per scenario, the drag sweep table, the
vortex line integrals and ``summary.json``. No plots or reports are generated.
"""

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from nsbinary.inspiral import Binary, build_reservoir, summarize_trace
from nsbinary.io import write_csv, write_json, write_trace_csv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config" / "inspiral.json")
    parser.add_argument("--background", type=Path, default=ROOT / "results")
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "inspiral")
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    background = args.background.resolve()
    for name in ("summary.json", "profiles.csv"):
        if not (background / name).exists():
            raise FileNotFoundError(
                f"Missing background file {background / name}; "
                "run scripts/compute_structure.py first.")
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    reservoir = build_reservoir(cfg, background)
    write_csv(out / "vortex_line_integrals.csv", reservoir.columns)
    scenarios = {}

    def run(name, this_cfg=cfg, rigid=False, initialization="locked",
            save=True, record=True):
        """Evolve one scenario, record its summary and optionally save the trace.

        Rebuild the line integrals when microphysics changes; sweep runs reuse
        the fiducial reservoir. Traces are decimated for export only.
        """
        res = reservoir if this_cfg is cfg else build_reservoir(this_cfg, background)
        trace = Binary(this_cfg, res, rigid=rigid).evolve(initialization)
        summary = {"rigid": rigid, "initialization": initialization,
                   **summarize_trace(trace)}
        if record:
            scenarios[name] = summary
        if save:
            write_trace_csv(out / f"{name}.csv", trace)
        print(f"{name}: heat={summary['heat_J']:.5e} J; "
              f"heat/orbit={summary['heat_over_orbital_binding_energy_loss']:.5e}",
              flush=True)
        return summary

    # Same original phase for both paper-start models: capture is an outcome.
    # Locked controls isolate later evolution conditional on being in resonance.
    run("fiducial_paper", initialization="paper")
    run("rigid_paper", rigid=True, initialization="paper")
    run("fiducial_locked")
    run("rigid_locked", rigid=True)
    run("permanently_pinned_locked", {**cfg, "unpinned_fraction": 0.})
    run("suppressed_mobility_locked", {**cfg, "unpinned_fraction": 1e-3})
    # Each sweep member has its own locked initialization, not the endpoint of
    # the previous trajectory. Fix the crust/core drag ratio, not the mean B.
    sweep = []
    for drag in cfg["drag_sweep_core_ratios"]:
        local = {**cfg, "core_drag_ratio": drag,
                 "inner_crust_drag_ratio": drag*cfg["drag_sweep_crust_over_core"]}
        summary = run(f"drag_{drag:g}", local, save=False, record=False)
        sweep.append({"core_drag_ratio": drag,
                      **{key: value for key, value in summary.items()
                         if key not in ("rigid", "initialization")}})
    write_csv(out / "drag_sweep.csv",
              {key: [row[key] for row in sweep] for key in sweep[0]})
    write_json(out / "summary.json",
               {"reservoir": reservoir.metadata, "scenarios": scenarios,
                "drag_sweep": sweep})
    print(f"Wrote results: {out}", flush=True)


if __name__ == "__main__":
    main()
