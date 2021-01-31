""" Functions for unit conversions in MMElemental """

__all__ = ["convert"]

from pint import UnitRegistry


def convert(quant, from_units, to_units):
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
    Any
        Converted quantity.
    """
    ureg = UnitRegistry()
    uquant = ureg.Quantity(quant, from_units)
    cquant, _ = uquant.to(to_units).to_tuple()
    return cquant
