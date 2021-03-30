from pydantic import Field, constr, validator
from mmelemental.models.forcefield.params import Params
from typing import Optional, Dict, Any, Union, List, Tuple
import qcelemental
import os
import pathlib

__all__ = ["Bonds"]


class Bonds(Params):
    lengths: qcelemental.models.types.Array[float] = Field(
        ..., description="Equilibrium bond lengths. Default unit is Angstroms."
    )
    lengths_units: Optional[str] = Field(
        "angstroms", description="Equilibrium bond lengths unit."
    )
    indices: List[Tuple[int, int, float]] = Field(  # type: ignore, need to make this field non-optional?
        ...,
        description="Particle indices for each bond and the bond order: (index1, index2, order).",
        min_items=1,
    )
    _path = os.path.join(pathlib.Path(__file__).parent.absolute(), "potentials", "*.py")

    def dict(self, *args, **kwargs):
        kwargs["exclude"] = {"provenance"}
        return super().dict(*args, **kwargs)
