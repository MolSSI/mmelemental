from pydantic import Field
from mmelemental.models.forcefield.params import Params
from mmelemental.util.units import LENGTH_DIM
from typing import Optional, List, Tuple, Union
from cmselemental.types import Array
import os
import pathlib

__all__ = ["Bonds"]


class Bonds(Params):
    lengths: Array[float] = Field(
        ..., description="Equilibrium bond lengths. Default unit is Angstroms."
    )
    lengths_units: Optional[str] = Field(
        "angstroms",
        description="Equilibrium bond lengths unit.",
        dimensionality=LENGTH_DIM,
    )
    connectivity: Optional[List[Tuple[Union[int, str], Union[int, str], float]]] = Field(  # type: ignore
        None,
        description="Particle indices  or names e.g. types for each bond and the bond order: (index1, index2, order).",
        min_items=1,
    )
    _path = os.path.join(pathlib.Path(__file__).parent.absolute(), "potentials", "*.py")
