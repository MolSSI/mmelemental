from pydantic import Field, validator
from typing import Optional
from mmelemental.models.base import ProtoModel
from cmselemental.types import Array

__all__ = ["Harmonic"]


class Harmonic(ProtoModel):
    """
    Linear spring angle model: Energy = 1/2 * spring * (angle - eq_angle)**2. "
    """

    spring: Optional[Array[float]] = Field(0, description="Angle spring constant. ")
    spring_units: Optional[str] = Field(
        "kJ/(mol*degrees**2)", description="Angle spring constant unit."
    )

    @validator("spring", allow_reuse=True)
    def _valid_length(cls, v):
        assert len(v.shape) == 1, "Angle spring constants must be a 1D array!"
        return v

    def dict(self, *args, **kwargs):
        kwargs["exclude"] = {"provenance"}
        return super().dict(*args, **kwargs)
