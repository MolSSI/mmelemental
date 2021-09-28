from pydantic import Field, validator, root_validator
from typing import Optional
from cmselemental.types import Array
from mmelemental.models.base import ProtoModel
from mmelemental.util.units import (
    LENGTH_DIM,
    MASS_DIM,
    TIME_DIM,
    SUBS_DIM,
    DIMENSIONLESS,
)

__all__ = ["Charmm"]


class Charmm(ProtoModel):
    """
    Charmm-style improper dihedral potential: Energy = energy * (1 + cos(periodicity * angle - phase)).
    """

    energy: Array[float] = Field(
        ...,
        description="Improper dihedral energy constant. Default unit is kJ/mol.",
    )
    energy_units: Optional[str] = Field(
        "kJ/mol",
        description="Improper dihedral energy constant unit.",
        dimensionality=MASS_DIM * LENGTH_DIM ** 2 / (TIME_DIM ** 2 * SUBS_DIM),
    )
    periodicity: Array[int] = Field(
        ...,
        description="Improper dihedral periodicity factor, must be >= 0.",
    )
    phase: Array[float] = Field(
        ...,
        description="Improper dihedral phase angle. Default unit is degrees.",
    )
    phase_units: Optional[str] = Field(
        "degrees",
        description="Improper dihedral phase angle unit.",
        dimensionality=DIMENSIONLESS,
    )

    @validator("energy", allow_reuse=True)
    def _valid_shape(cls, v):
        assert (
            len(v.shape) == 1
        ), "Improper dihedral energy constant must be a 1D array!"
        return v

    @validator("periodicity", allow_reuse=True)
    def _valid_periodicity(cls, v):
        assert (v >= 0).all(), "Improper dihedral periodicity must be >= 0."
        return v

    @root_validator(allow_reuse=True)
    def _valid_arrays(cls, values):
        energy_shape = values["energy"].shape
        periodicity_shape = values["periodicity"].shape
        phase_shape = values["phase"].shape
        assert (
            energy_shape == periodicity_shape == phase_shape
        ), f"Energy ({energy_shape}), periodocity ({periodicity_shape}), and phase ({phase_shape}) must have the same shape."
        return values
