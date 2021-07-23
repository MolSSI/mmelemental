from pydantic import Field
from mmelemental.models.forcefield.params import Params
from typing import Optional, List, Tuple, Union
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
    connectivity: List[Tuple[Union[int, str], Union[int, str], float]] = Field(  # type: ignore
        ...,
        description="Particle indices  or names e.g. types for each bond and the bond order: (index1, index2, order).",
        min_items=1,
    )
    _path = os.path.join(pathlib.Path(__file__).parent.absolute(), "potentials", "*.py")

    def dict(self, *args, **kwargs):
        kwargs["exclude"] = {"provenance"}
        return super().dict(*args, **kwargs)
