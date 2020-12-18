from .base import SimInput
from mmelemental.models.molecule.mm_molecule import Molecule
from pydantic import Field
from typing import Any, Tuple, Union

class OptimInput(SimInput):
    mol: Union[Tuple[Molecule], Molecule] = Field(..., description = 'Molecular mechanics molecule object(s).')
    tol: float = Field(None, description = 'Tolerance used to indicate when the optimization scheme has converd.')