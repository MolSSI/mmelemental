from .base import SimBase
from mmelemental.models.molecule.mm_molecule import Molecule
from pydantic import Field
from typing import Any, Tuple, Union

class Optimization(SimBase):
    tol: float = Field(None, description = 'Tolerance used to indicate when the optimization scheme has converd.')