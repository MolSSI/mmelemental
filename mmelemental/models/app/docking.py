from mmelemental.models.app.base import SimInput
from mmelemental.models.base import Base
from mmelemental.models.molecule import Molecule
from mmelemental.models.collect import Ensemble, Trajectory
from pydantic import Field
from typing import List, Dict, Optional, Tuple, Union


class DockInput(SimInput):
    ligand: Union[List[Molecule], Molecule] = Field(
        ...,
        description="Ligand molecule object such as a drug. See the :class:``Molecule`` class. ",
    )
    receptor: Union[List[Molecule], Molecule] = Field(
        ...,
        description="Receptor molecule object such as a protein. See the :class:``Molecule`` class. ",
    )
    searchSpace: Optional[Tuple[float, float, float, float, float, float]] = Field(
        None,
        description="A 3D box defined by (xmin, xmax, ymin, ymax, zmin, zmax)."
        "The search space effectively restricts where the movable atoms, including those in the flexible side chains, should lie.",
    )


class DockOutput(Base):
    dockInput: DockInput = Field(..., description="Docking input model.")
    ligands: List[Molecule] = Field(
        ...,                                                                                                                                
        description="Molecule object(s) representing ligands. See the :class:``Molecule`` class.",
    )
    receptors: Optional[List[Molecule]] = Field(
        None,
        description="Molecule object(s) representing flexible receptors e.g. side-chains. See the :class:``Molecule`` class.",
    )
    scores: List[float] = Field(
        ...,
        description="Scores for each pose. Typically the score is the binding affinity. Length must be equal to the number of poses.",
    )
    scores_units: Optional[str] = Field(
        "kJ/mol",
        description="Score function unit."
    )
    observables: Optional[Dict[str, List[float]]] = Field(
        None,
        description="Physical observables for each pose. E.g. observables={'RMSD':[...]}."
    )
    observables_units: Optional[Dict[str,str]] = Field(
        None,
        description="Physical observables units. E.g. observables_units={'RMSD':'angstrom'}."
    )
