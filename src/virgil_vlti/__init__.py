"""VLTI calibration for virgil: GRAVITY, plus MATISSE headers and a PIONIER pair test."""

from importlib.metadata import PackageNotFoundError, version

from .calibrate import calibrate_exposure, flux_ratio, frame_times, group_by_night
from .load import to_oidata
from .products import Exposure, classify, read_exposure
from .write import write_calibrated

try:
    __version__ = version("virgil-vlti")
except PackageNotFoundError:  # running from a source tree that is not installed
    __version__ = "unknown"

__all__ = [
    "__version__",
    "Exposure",
    "calibrate_exposure",
    "classify",
    "flux_ratio",
    "frame_times",
    "group_by_night",
    "read_exposure",
    "to_oidata",
    "write_calibrated",
]
