from pydantic import Field, constr, validator
from mmelemental.models.forcefield.params import Params
import qcelemental
import hashlib
from typing import Optional, Dict, Any, Union, List
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

    # Validators
    @validator("angles")
    def _angles_length(cls, v, values):
        assert len(v.shape) == 1, "Bond lengths must be a 1D array!"
        return v
