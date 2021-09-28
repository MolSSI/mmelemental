from pydantic import Field, validator, root_validator
from typing import Optional, List
from cmselemental.types import Array
from mmelemental.models.base import ProtoModel
from mmelemental.util.units import (
    LENGTH_DIM,
    MASS_DIM,
    TIME_DIM,
    SUBS_DIM,
    DIMENSIONLESS,
)

__all__ = ["CharmmMulti"]


class CharmmMulti(ProtoModel):
    """
    Charmm-style multiple improper dihedral potential terms: Energy = SUM[ energy * (1 + cos(periodicity * angle - phase)) ].
    """

    energy: List[Array[float]] = Field(
        ...,
        description="Improper dihedral energy constant. Default unit is kJ/mol.",
    )
    energy_units: Optional[str] = Field(
        "kJ/mol",
        description="Improper dihedral energy constant unit.",
        dimensionality=MASS_DIM * LENGTH_DIM ** 2 / (TIME_DIM ** 2 * SUBS_DIM),
    )
    periodicity: List[Array[int]] = Field(
        ...,
        description="Improper dihedral periodicity factor, must be >= 0.",
    )
    phase: List[Array[float]] = Field(
        ...,
        description="Improper dihedral phase angle. Default unit is degrees.",
    )
    phase_units: Optional[str] = Field(
        "degrees",
        description="Improper dihedral phase angle unit.",
        dimensionality=DIMENSIONLESS,
    )

    @validator("energy", "periodicity", "phase", allow_reuse=True)
    def _valid_shape(cls, v):
        assert v[-1].shape, "Array must be 2D!"
        return v

    @validator("periodicity", allow_reuse=True)
    def _valid_periodicity(cls, v):
        for arr in v:
            assert (arr >= 0).all(), "Improper dihedral periodicity must be >= 0."
        return v

    @root_validator(allow_reuse=True)
    def _valid_arrays(cls, values):
        energy_len = len(values["energy"])
        periodicity_len = len(values["periodicity"])
        phase_len = len(values["phase"])
        assert (
            energy_len == periodicity_len == phase_len
        ), f"Energy ({energy_len}), periodocity ({periodicity_len}), and phase ({phase_len}) must have the same shape."
        return values
