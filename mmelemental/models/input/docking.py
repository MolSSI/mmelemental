from .base import SimBase
from typing import List, Optional, Tuple, Union
from mmelemental.models.molecule.mm_molecule import Molecule
from pydantic import Field

class DockingInput(SimBase):
    ligand: Molecule = Field(
        ..., description="Molecule model for a candidate ligand (e.g. drug)."
    )
    receptor: Molecule = Field(
        ..., description="Molecule model for a receptor (e.g. protein)."
    )
    searchSpace: Optional[List[Tuple[float, float, float, float, float, float]]] = Field(
        None,
        description="A 3D box defined by (xmin, xmax, ymin, ymax, zmin, zmax)."
        "The search space effectively restricts where the movable atoms, including those in the flexible side chains, should lie.",
    )