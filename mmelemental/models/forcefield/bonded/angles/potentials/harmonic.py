from pydantic import Field, root_validator
from typing import Optional
import qcelemental
from ..angle_params import AngleParams

__all__ = ["Harmonic"]


class Harmonic(AngleParams):
    """ 
    Linear spring angle model: Energy = 1/2 * spring * (angle - eq_angle)**2. "
    """

    spring: Optional[qcelemental.models.types.Array[float]] = Field(
        0, description="Angle spring constant. "
    )
    spring_units: Optional[str] = Field(
        "kJ/(mol*degrees**2)", description="Angle spring constant unit."
    )

    @root_validator
    def _valid_length(cls, values):
        assert (
            len(values["spring"].shape) == 1
        ), "Angle spring constants must be a 1D array!"
        return values
