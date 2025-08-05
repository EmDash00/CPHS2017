from enum import IntEnum
from math import tau

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import circmean


def wrap_angle(
    x: np.ndarray[tuple[int, ...], np.dtype[np.float64]],
    out: np.ndarray[tuple[int, ...], np.dtype[np.float64]] | None = None,
    **kwargs,
):
    return np.mod(x, tau, out, **kwargs)


def angular_clip(
    x: np.ndarray[tuple[int, ...], np.dtype[np.float64]],
    lower: float = 0,
    higher: float = tau,
    out: np.ndarray[tuple[int, ...], np.dtype[np.float64]] | None = None,
    **kwargs,
):
    out = out if out is not None else np.empty_like(x)
    wrap_angle(x - lower, out=out)
    out *= (higher - lower) / tau
    out -= lower

    return out

def angular_diff(a, b):
    return (a - b + np.pi) % tau - np.pi


def circmedian(angles):
    """Fast approximate circular median via mean + local optimization."""
    # Step 1: Compute circular mean as initial guess

    mean_angle = circmean(angular_clip(angles))

    # Step 2: Minimize sum of circular deviations near the mean
    def sum_circular_deviations(theta):
        return np.sum(np.abs(angular_diff(theta, angles)))

    return minimize_scalar(
        sum_circular_deviations,
        bounds=(mean_angle - 0.25 * tau, mean_angle + 0.25 * tau),
        method="bounded",
    ).x


class StatArray(np.ndarray):

    def __new__(cls, input_array):
        return np.asarray(input_array).view(cls)

    def __array_finalize__(self, obj):
        if obj is None: 
            return
        
    @property
    def x(self):
        return self[0]


PATTERN = (
    r"(?P<date>\d+-\d+)_"
    r"(?P<protocol_name>[a-z0-9]+)_"
    r"(?P<order>so|fo)_"
    r"r-(?P<ref_str>[a-z-+AEO]+)_"
    r"d-(?P<dist_str>[a-z-+AEO]+)_"
    r"s(?P<seed>\d+)_"
    r"p(?P<shift_id>\d)"
    r"a(?P<assay_index>\d)_"
    r"(?P<sign>[+-]1)"
)


class SignalType(IntEnum):
    ZERO = 0
    SUM_OF_SINES_PLUS_EVEN = 1
    SUM_OF_SINES_PLUS_ODD = 2
    SUM_OF_SINES_PLUS_ALL = 3
    SUM_OF_SINES_MINUS_EVEN = -1
    SUM_OF_SINES_MINUS_ODD = -2
    SUM_OF_SINES_MINUS_ALL = -3


SIGNAL_TYPE_STR_MAP = {
    "zer": SignalType.ZERO,
    "sos+E": SignalType.SUM_OF_SINES_PLUS_EVEN,
    "sos+O": SignalType.SUM_OF_SINES_PLUS_ODD,
    "sos+A": SignalType.SUM_OF_SINES_PLUS_ALL,
    "sos-E": SignalType.SUM_OF_SINES_MINUS_EVEN,
    "sos-O": SignalType.SUM_OF_SINES_MINUS_ODD,
    "sos-A": SignalType.SUM_OF_SINES_MINUS_ALL,
}


class AssayType(IntEnum):
    """
    Codes the assay type as a bit vector using 6 bits.
    Reference information is in 3 most significant bits.
    Disturbance information is in 3 least significant bits.

    In each group of 3 bits, the following information is encoded:

    Bit 1 (MSB): Sign bit (1 is negative and 0 is positive)
    Bit 2: Has even frequencies
    Bit 3: Has odd frequencies
    """

    SUM_OF_SINES_R_PLUS_ALL_D_ZERO = 0b011_000
    SUM_OF_SINES_R_ZERO_D_PLUS_ALL = 0b000_111
    SUM_OF_SINES_R_PLUS_EVEN_D_PLUS_ODD = 0b010_001
    SUM_OF_SINES_R_PLUS_ODD_D_PLUS_EVEN = 0b001_010
    SUM_OF_SINES_R_MINUS_EVEN_D_PLUS_ODD = 0b110_001
    SUM_OF_SINES_R_MINUS_ODD_D_PLUS_EVEN = 0b101_010


ASSAY_TYPE_ARR = list(AssayType)


def assay_index_to_type(index: int) -> AssayType:
    if index < len(AssayType):
        return ASSAY_TYPE_ARR[index]

    if index == 6:
        return AssayType.SUM_OF_SINES_R_PLUS_EVEN_D_PLUS_ODD

    raise ValueError("Invalid assay index.")


ORDER_STRING_TO_NUM = {"fo": 1, "so": 2}
