from mmelemental.models.forcefield.params import Params
from pydantic import Field
from typing import Optional
import os
import pathlib

__all__ = ["NonBonded"]


class NonBonded(Params):
    """Model that describes non-bonded interactions between particles."""

    _path = os.path.join(pathlib.Path(__file__).parent.absolute(), "potentials", "*.py")
    combination_rule: Optional[str] = Field(
        "Lorentz-Berthelot", description="Combination rule for the force field."
    )
