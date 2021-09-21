from pydantic import Field, validator
from typing import Optional
from mmelemental.models.base import ProtoModel
from mmelemental.util.units import (
    MASS_DIM,
    TIME_DIM,
    SUBS_DIM,
)
from cmselemental.types import Array

__all__ = ["Gromos96"]


class Gromos96(ProtoModel):
    """
    GROMOS-96 spring bond model: Energy = 1/4 * spring * (distance**2 - length**2)**2. "
    """

    spring: Array[float] = Field(
        ..., description="Bond spring constant. Default unit is kJ/(mol*angstrom**2)."
    )
    spring_units: Optional[str] = Field(
        "kJ/(mol*angstrom**2)",
        description="Bond spring constant unit.",
        dimensionality=MASS_DIM / (TIME_DIM ** 2 * SUBS_DIM),
    )

    @validator("spring", allow_reuse=True)
    def _valid_length(cls, v):
        assert len(v.shape) == 1, "Bond spring constants must be a 1D array!"
        return v
