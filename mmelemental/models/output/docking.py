from typing import List, Optional
from mmelemental.models.input.docking import DockingInput
from mmelemental.models.output.sim import SimOutput
from mmelemental.models.base import Base
from pydantic import Field


class DockingOutput(Base):
    dockingInput: DockingInput = Field(..., description="Docking input model.")
    ligand: SimOutput = Field(
        ...,
        description="Simulation output for the ligand, including its pose and score.",
    )
    receptor: Optional[SimOutput] = Field(
        None,
        description="Simulation output for non-rigid receptors i.e. conformation and orientation of the flexible side chains in the receptor relative to the ligand.",
    )


class AffinityOutput(Base):
    Docking_Output: DockingOutput
    Affinity: float
