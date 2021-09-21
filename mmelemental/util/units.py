""" Functions for unit conversions in MMElemental """

__all__ = ["convert", "Quantity"]

from typing import Any
from pint import UnitRegistry, Quantity, util

LENGTH_DIM = util.UnitsContainer({"[length]": 1})
MASS_DIM = util.UnitsContainer({"[mass]": 1})
TIME_DIM = util.UnitsContainer({"[time]": 1})
CURRENT_DIM = util.UnitsContainer({"[current]": 1})
SUBS_DIM = util.UnitsContainer({"[substance]": 1})
DIMENSIONLESS = util.UnitsContainer({})


def convert(quant: Any, from_units: str, to_units: str):
    """Converts quantity from units in 'from_units' to units in 'to_units'.
    All units supported by pint are supoorted in MMSchema.
    Parameters
    ----------
    quant: Any
        Quantity to convert.
    from_units: str
        Units to convert from.
    to_units: str
        Units to convert to.
    Returns
    --------
    numpy.ndarray, or float, or int
        Converted quantity.
    """
    if quant is not None:
        ureg = UnitRegistry()
        uquant = ureg.Quantity(quant, from_units)
        cquant, _ = uquant.to(to_units).to_tuple()
    else:
        cquant = None
    return cquant
