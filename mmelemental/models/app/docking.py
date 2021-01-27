from mmelemental.models.app.base import SimInput
from mmelemental.models.base import Base
from mmelemental.models.molecule.mm_mol import Mol
from mmelemental.models.trajectory.mm_traj import Traj
from pydantic import Field
from typing import List, Optional, Tuple


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


class DockingOutput(Base):
    dockingInput: DockingInput = Field(..., description="Docking input model.")
    ligand: Traj = Field(
        ...,
        description="Simulation output for the ligand, including its pose and score.",
    )
    receptor: Optional[Traj] = Field(
        None,
        description="Simulation output for non-rigid receptors i.e. conformation and orientation of the flexible side chains in the receptor relative to the ligand.",
    )


class AffinityOutput(Base):
    Docking_Output: DockingOutput
    Affinity: float
