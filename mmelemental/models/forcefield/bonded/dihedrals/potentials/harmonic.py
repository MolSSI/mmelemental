from pydantic import Field, validator
from typing import Optional
from cmselemental.types import Array
from mmelemental.models.base import ProtoModel
from mmelemental.util.units import (
    LENGTH_DIM,
    MASS_DIM,
    TIME_DIM,
    SUBS_DIM,
)

__all__ = ["Harmonic"]


class Harmonic(ProtoModel):
    """
    Simple periodic dihedral potential: Energy = energy * (1 + sign * cos(periodicity * angle)).
    """

    energy: Array[float] = Field(
        ...,
        description="Dihedral energy constant. Default unit is kJ/mol.",
    )
    energy_units: Optional[str] = Field(
        "kJ/mol",
        description="Dihedral energy constant unit.",
        dimensionality=MASS_DIM * LENGTH_DIM ** 2 / (TIME_DIM ** 2 * SUBS_DIM),
    )
    periodicity: Array[int] = Field(
        ...,
        description="Dihedral periodicity term, must be >= 0.",
    )
    sign: Optional[int] = Field(
        1,
        description="Multiplication factor for the cosine term. Must be either -1 or +1.",
    )

    @validator("energy", allow_reuse=True)
    def _valid_shape(cls, v):
        assert len(v.shape) == 1, "Dihedral energy constant must be a 1D array!"
        return v

    @validator("sign", allow_reuse=True)
    def _valid_sign(cls, v):
        assert v == 1 or v == -1, "Dihedral sign term can be either +1 or -1."
        return

    @validator("periodicity", allow_reuse=True)
    def _valid_periodicity(cls, v):
        assert (v >= 0).all(), "Dihedral periodicity must be >= 0."
        return v
