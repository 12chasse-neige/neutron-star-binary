"""Writers for the numerical result files.

Only aligned CSV columns and JSON metadata are produced; the project contains no
plotting or report generation.
"""

import csv
import json
from pathlib import Path

import numpy as np


def write_csv(path, columns):
    """Write dictionary columns in insertion order over the full supplied grid."""
    with Path(path).open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(zip(*columns.values()))


def write_trace_csv(path, columns, max_rows=15000):
    """Write a possibly decimated trace while spanning the full sampled range.

    Retained rows are not uniformly spaced in time or frequency because the
    adaptive solver nodes are included before decimation.
    """
    count = len(next(iter(columns.values())))
    rows = (np.unique(np.linspace(0, count - 1, min(count, max_rows), dtype=int))
            if max_rows else range(count))
    with Path(path).open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(zip(*(np.asarray(v)[rows] for v in columns.values())))


def write_json(path, obj):
    """Write JSON, converting NumPy scalars and rejecting NaN/Inf."""
    def convert(value):
        return value.item() if isinstance(value, np.generic) else str(value)
    Path(path).write_text(json.dumps(obj, indent=2, default=convert,
                                     allow_nan=False) + "\n")
