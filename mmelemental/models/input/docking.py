from .base import SimInput
from typing import List, Optional, Tuple, Union
from ..molecule.mm_mol import Mol
from pydantic import Field


class DockingInput(SimInput):
    ligand: Mol = Field(
        ...,
        description="Molecule model for a candidate ligand (e.g. drug). See  the :class:``Mol``.",
    )
    receptor: Mol = Field(
        ...,
        description="Molecule model for a receptor (e.g. protein). See  the :class:``Mol``.",
    )
    searchSpace: Optional[
        List[Tuple[float, float, float, float, float, float]]
    ] = Field(
        None,
        description="A 3D box defined by (xmin, xmax, ymin, ymax, zmin, zmax)."
        "The search space effectively restricts where the movable atoms, including those in the flexible side chains, should lie.",
    )
