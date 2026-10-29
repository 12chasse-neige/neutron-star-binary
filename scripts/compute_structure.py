#!/usr/bin/env python3
"""Solve the APR neutron-star background and write numeric structure results.

Run with ``.venv/bin/python scripts/compute_structure.py``; default paths resolve
from the repository root. Outputs in ``--output`` (default ``results/``):
``summary.json``, ``profiles.csv`` and ``layers.csv``. No plots or reports are
generated.
"""

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from nsbinary.eos import load_apr
from nsbinary.io import write_csv, write_json
from nsbinary.structure import build_summary, solve_target


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config" / "structure.json")
    parser.add_argument("--data", type=Path, default=ROOT / "data")
    parser.add_argument("--output", type=Path, default=ROOT / "results")
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    if not 0 < cfg["neutron_drip_rest_mass_density_g_cm3"] < cfg["core_crust_rest_mass_density_g_cm3"]:
        raise ValueError("Require positive neutron-drip density below the core-crust density.")
    eos = load_apr(args.data)
    star = solve_target(eos, cfg["target_mass_msun"],
                        surface_pseudo_enthalpy=cfg["surface_pseudo_enthalpy"],
                        rtol=cfg["relative_tolerance"],
                        max_step_km=cfg["max_step_km"])
    summary = build_summary(star, cfg)
    profile = star.profile(cfg["profile_points"],
                           [summary[key]["radius_km"] for key in
                            ("core_crust_boundary", "neutron_drip_boundary")])
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    write_json(out / "summary.json", summary)
    write_csv(out / "profiles.csv", profile)
    write_csv(out / "layers.csv",
              {key: [layer[key] for layer in summary["layers"]]
               for key in summary["layers"][0]})
    print(f"Solved: M={star.mass_msun:.9f} Msun, R={star.radius_km:.6f} km -> {out}",
          flush=True)


if __name__ == "__main__":
    main()
