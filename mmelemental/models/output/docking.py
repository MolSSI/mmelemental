from typing import List, Optional
from mmelemental.models.input.docking import DockingInput
from mmelemental.models.molecule.mm_molecule import Molecule
from mmelemental.models.base import Base
from pydantic import Field

class DockingOutput(Base):
    dockingInput: DockingInput = Field(..., description="Docking input model.")
    poses: List[Molecule] = Field(
        ...,
        description="Conformation and orientation of the candidate ligand relative to the receptor.",
    )
    scores: List[float] = Field(
        ...,
        description="A metric for evaluating a particular pose. Length of scores must be equal to length of poses.",
    )
    flexible: Optional[List[Molecule]] = Field(
        None,
        description="Conformation and orientation of the flexible side chains in the receptor relative to the ligand.",
    )

class AffinityOutput(Base):
    Docking_Output: DockingOutput
    Affinity: float