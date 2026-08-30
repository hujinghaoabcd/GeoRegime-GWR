"""Research implementation of Geo-Regime GWR.

The public surface is intentionally small because this repository is for
method development and paper experiments rather than a mature software API.
"""

from .gwr import BasicGWR, GWRResult
from .grgwr import GRGWRBaseline, GRGWRResult

__all__ = ["BasicGWR", "GWRResult", "GRGWRBaseline", "GRGWRResult"]
