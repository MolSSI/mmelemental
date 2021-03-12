from mmelemental.models.forcefield.params import Params
from pydantic import Field, validator
from typing import Optional
import os
import pathlib

__all__ = ["NonBondedParams"]


class NonBondedParams(Params):

    name: Optional[str] = Field(None, description="Name or form of the potential.")
    _path_name = os.path.join(
        pathlib.Path(__file__).parent.absolute(), "potentials", "*.py"
    )

    # Validators
    @validator("name", always=True)
    def _set_name(cls, v, values):
        if v is not None:
            assert v == cls.__name__
        return cls.__name__
