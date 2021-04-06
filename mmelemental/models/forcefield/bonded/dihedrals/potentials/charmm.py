from pydantic import Field, validator
from typing import Optional
import qcelemental
from mmelemental.models.base import ProtoModel

__all__ = ["Charmm"]


class Charmm(ProtoModel):
    """
    Charmm-style dihedral potential: Energy = energy * (1 + cos(periodicity * angle - phase)).
    """

    energy: qcelemental.models.types.Array[float] = Field(
        ...,
        description="Dihedral energy constant. Default unit is kJ/mol.",
    )
    energy_units: Optional[str] = Field(
        "kJ/mol", description="Dihedral energy constant unit."
    )
    periodicity: qcelemental.models.types.Array[int] = Field(
        ...,
        description="Dihedral periodicity term, must be >= 0.",
    )
    phase: qcelemental.models.types.Array[float] = Field(
        ...,
        description="Dihedral phase angle. Default unit is degrees.",
    )
    phase_units: Optional[str] = Field(
        "degrees", description="Dihedral phase angle unit."
    )

    @validator("energy", allow_reuse=True)
    def _valid_shape(cls, v):
        assert len(v.shape) == 1, "Dihedral energy constant must be a 1D array!"
        return v

    @validator("periodicity", allow_reuse=True, each_item=True)
    def _valid_periodicity(cls, v):
        assert v >= 0, "Dihedral periodicity must be >= 0."
        return

    def dict(self, *args, **kwargs):
        kwargs["exclude"] = {"provenance"}
        return super().dict(*args, **kwargs)
