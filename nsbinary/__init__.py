"""Neutron-star structure and two-component inspiral core models.

Submodules
----------
units      Physical constants and SI/geometrized conversions.
eos        Cold tabulated barotrope adapter (APR) and monotone resampling.
structure  TOV, Hartle slow rotation, static tides and layer classification.
inspiral   Reservoir line integrals and the two-rotator inspiral model.
io         CSV/JSON writers for numerical results.
"""

from . import eos, inspiral, io, structure, units

__all__ = ["eos", "inspiral", "io", "structure", "units"]
