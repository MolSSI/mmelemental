from pydantic import Field, validator
from mmelemental.models.forcefield.params import Params
from cmselemental.types import Array
from typing import Optional, List, Tuple
import os
import pathlib

__all__ = ["Dihedrals"]


class Dihedrals(Params):
    _path = os.path.join(pathlib.Path(__file__).parent.absolute(), "potentials", "*.py")
    angles: Array[float] = Field(
        None, description="Equilibrium dihedral angles. Default unit is degrees."
    )
    angles_units: Optional[str] = Field(
        "degrees", description="Equilibrium dihedral angle units."
    )
    connectivity: List[Tuple[int, int, int, int]] = Field(  # type: ignore
        ...,
        description="Particle indices for each dihedral angle.",
        min_items=1,
    )
    weights: Optional[Array[float]] = Field(
        None,
        description="Something to consider later on? Ses CHARMM dihedral_style for LAMMPS.",
    )

    # Validators
    @validator("angles")
    def _angles_shape(cls, v, values):
        assert len(v.shape) == 1, "Bond lengths must be a 1D array!"
        return v

    @validator("weights")
    def _valid_weights(cls, v):
        import numpy

        unique_weights = numpy.unique(v)
        assert numpy.all(
            [weight in (0, 0.5, 1) for weight in unique_weights]
        ), "Weight factor can be 0, 0.5, or 1."
        return v
