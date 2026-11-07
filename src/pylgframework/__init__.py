"""pylgframework: Local Group barycentre with propagated observational uncertainty."""
from . import fiducials
from .barycentre import Result, barycentre
from .compare import Comparison, compare
from .sampling import DEFAULT_SIGMAS

__version__ = "0.1.0"
__all__ = ["barycentre", "compare", "fiducials", "Result", "Comparison", "DEFAULT_SIGMAS"]
