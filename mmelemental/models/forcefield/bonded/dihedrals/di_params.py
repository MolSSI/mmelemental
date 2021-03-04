from mmelemental.models.forcefield.params import Params
import qcelemental
from pydantic import Field, validator
from typing import Optional
import os
import pathlib

__all__ = ["DihedralParams"]


class DihedralParams(Params):
    angles: qcelemental.models.types.Array[float] = Field(
        None, description="Equilibrium dihedral angles. Default unit is degrees."
    )
    angles_units: Optional[str] = Field(
        "degrees", description="Equilibrium dihedral angle units."
    )
    name: Optional[str] = Field(None, description="Name or form of the potential.")
    _path_name = os.path.join(
        pathlib.Path(__file__).parent.absolute(), "potentials", "*.py"
    )

    # Validators
    @validator("angles")
    def _angles_length(cls, v, values):
        assert len(v.shape) == 1, "Bond lengths must be a 1D array!"
        return v

    @validator("name", always=True)
    def _set_name(cls, v, values):
        if v is not None:
            assert v == cls.__name__
        return cls.__name__
