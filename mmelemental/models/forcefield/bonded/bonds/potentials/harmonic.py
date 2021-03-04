from pydantic import Field, root_validator
from typing import Optional
import qcelemental
from ..bond_params import BondParams

__all__ = ["Harmonic"]


class Harmonic(BondParams):
    """
    Linear spring bond model: Energy = 1/2 * spring * (distance - length)**2. "
    """

    spring: qcelemental.models.types.Array[float] = Field(
        0, description="Bond spring constant. Default unit is kJ/(mol*angstrom**2)."
    )
    spring_units: Optional[str] = Field(
        "kJ/(mol*angstrom**2)", description="Bond spring constant unit."
    )

    @root_validator
    def _valid_length(cls, values):
        assert (
            len(values["spring"].shape) == 1
        ), "Bond spring constants must be a 1D array!"
        return values
