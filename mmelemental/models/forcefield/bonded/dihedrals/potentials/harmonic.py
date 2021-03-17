from pydantic import Field, validator
from typing import Optional
import qcelemental
from mmelemental.models.base import ProtoModel

__all__ = ["Harmonic"]


class Harmonic(ProtoModel):
    """
    Linear spring dihedral model: Energy = 1/2 * spring * (angle - eq_angle)**2. "
    """

    spring: qcelemental.models.types.Array[float] = Field(
        ...,
        description="Dihedral spring constant. Default unit is kJ/(mol*degrees**2).",
    )
    spring_units: Optional[str] = Field(
        "kJ/(mol*degrees**2)", description="Dihedral spring constant unit."
    )

    @validator("spring", allow_reuse=True)
    def _valid_length(cls, v):
        assert len(v.shape) == 1, "Bond spring constants must be a 1D array!"
        return v

    def dict(self, *args, **kwargs):
        kwargs["exclude"] = {"provenance"}
        return super().dict(*args, **kwargs)
