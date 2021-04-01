from pydantic import Field, validator
from mmelemental.models.forcefield.params import Params
import qcelemental
from typing import Optional, List, Tuple
import os
import pathlib

__all__ = ["Dihedrals"]


class Dihedrals(Params):
    _path = os.path.join(pathlib.Path(__file__).parent.absolute(), "potentials", "*.py")
    angles: qcelemental.models.types.Array[float] = Field(
        None, description="Equilibrium dihedral angles. Default unit is degrees."
    )
    angles_units: Optional[str] = Field(
        "degrees", description="Equilibrium dihedral angle units."
    )
    indices: List[Tuple[int, int, int, int]] = Field(  # type: ignore, need to make this field non-optional?
        ...,
        description="Particle indices for each dihedral angle.",
        min_items=1,
    )

    # Validators
    @validator("angles")
    def _angles_length(cls, v, values):
        assert len(v.shape) == 1, "Bond lengths must be a 1D array!"
        return v
