from pydantic import Field, root_validator
from typing import Optional
import qcelemental
from ..di_params import DihedralsParams

__all__ = ["Harmonic"]


class Harmonic(DihedralsParams):
    """ 
    Linear spring dihedral model: Energy = 1/2 * spring * (angle - eq_angle)**2. "
    """

    spring: qcelemental.models.types.Array[float] = Field(
        ..., description="Dihedral spring constant. Default unit is kJ/(mol*degrees**2)."
    )
    spring_units: Optional[str] = Field(
        "kJ/(mol*degrees**2)", description="Dihedral spring constant unit."
    )

    @root_validator
    def _valid_length(cls, values):
        assert (
            len(values["spring"].shape) == 1
        ), "Bond spring constants must be a 1D array!"
        return values
