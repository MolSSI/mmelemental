from mmelemental.models.app.base import SimInput
from mmelemental.models.base import Base
from mmelemental.models.molecule import Molecule
from mmelemental.models.collect import Ensemble, Trajectory
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
    dockInput: DockInput = Field(..., description="Docking input model.")
    ensemble: Ensemble = Field(
        ...,
        description="Ensemble output for the ligand pose and score, and optionally the (flexible) receptor.",
    )


class AffinityOutput(Base):
    dockOutput: DockOutput
    affinity: float
