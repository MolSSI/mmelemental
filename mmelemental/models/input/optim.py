from .base import SimInput
from ..molecule.mm_mol import Mol
from pydantic import Field
from typing import Any, Tuple, Union


class OptimInput(SimInput):
    mol: Union[Tuple[Mol], Mol] = Field(
        ...,
        description="Molecular mechanics molecule object(s). See  the :class:``Mol``.",
    )
    tol: float = Field(
        None,
        description="Tolerance used to indicate when the optimization scheme has converd.",
    )
