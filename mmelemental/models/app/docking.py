from mmelemental.models.app.base import SimInput
from mmelemental.models.base import Base
from mmelemental.models.molecule.mm_mol import Molecule
from mmelemental.models.trajectory.mm_traj import Trajectory
from pydantic import Field
from typing import List, Optional, Tuple


class DockInput(SimInput):
    searchSpace: Optional[
        List[Tuple[float, float, float, float, float, float]]
    ] = Field(
        None,
        description="A 3D box defined by (xmin, xmax, ymin, ymax, zmin, zmax)."
        "The search space effectively restricts where the movable atoms, including those in the flexible side chains, should lie.",
    )


class DockOutput(Base):
    dockingInput: DockInput = Field(..., description="Docking input model.")
    ligand: Trajectory = Field(
        ...,
        description="Simulation output for the ligand, including its pose and score.",
    )
    receptor: Optional[Trajectory] = Field(
        None,
        description="Simulation output for non-rigid receptors i.e. conformation and orientation of the flexible side chains in the receptor relative to the ligand.",
    )


class AffinityOutput(Base):
    dockOutput: DockOutput
    affinity: float
