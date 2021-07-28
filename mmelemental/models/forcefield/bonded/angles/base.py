from pydantic import Field, validator
from mmelemental.models.forcefield.params import Params
from cmselemental.types import Array
from typing import Optional, Tuple, List, Union
import os
import pathlib

__all__ = ["Angles"]


class Angles(Params):
    _path = os.path.join(pathlib.Path(__file__).parent.absolute(), "potentials", "*.py")
    angles: Array[float] = Field(
        ..., description="Equilibrium angles. Default unit is degrees."
    )
    angles_units: Optional[str] = Field(
        "degrees", description="Equilibrium angle units."
    )
    connectivity: List[Tuple[Union[int, str], Union[int, str], Union[int, str]]] = Field(  # type: ignore, need to make this field non-optional?
        ...,
        description="Particle indices for each angle.",
        min_items=1,
    )

    # Validators
    @validator("angles")
    def _angles_length(cls, v, values):
        assert len(v.shape) == 1, "Bond lengths must be a 1D array!"
        return v

    def dict(self, *args, **kwargs):
        kwargs["exclude"] = {"provenance"}
        return super().dict(*args, **kwargs)
