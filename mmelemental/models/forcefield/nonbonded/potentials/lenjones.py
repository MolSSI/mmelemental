from pydantic import Field, root_validator
from typing import Optional
from mmelemental.models.base import ProtoModel
from mmelemental.util.units import (
    LENGTH_DIM,
    MASS_DIM,
    TIME_DIM,
    SUBS_DIM,
)
from cmselemental.types import Array

__all__ = ["LennardJones"]


class LennardJones(ProtoModel):

    epsilon: Array[float] = Field(
        ...,
        description="The epsilon (well depth) Lennard-Jones parameter. Default unit is kJ/mol.",
    )
    epsilon_units: Optional[str] = Field(
        "kJ/mol",
        description="Units for the Lennard-Jones epsilon (well depth) constant.",
        dimensionality=MASS_DIM * LENGTH_DIM ** 2 / (TIME_DIM ** 2 * SUBS_DIM),
    )
    sigma: Array[float] = Field(
        ...,
        description="The distance at which the Lennard-Jones potential is 0. Default unit is angstroms.",
    )
    sigma_units: Optional[str] = Field(
        "angstrom",
        description="Units for the Lennard-Jones sigma constant.",
        dimensionality=LENGTH_DIM,
    )

    @root_validator(allow_reuse=True)
    def _valid_length(cls, values):
        assert len(values["epsilon"].shape) == 1, "epsilon must be a 1D array!"
        assert len(values["sigma"].shape) == 1, "sigma must be a 1D array!"
        assert len(values["epsilon"]) == len(
            values["sigma"]
        ), "epsilon and sigma must be of equal length!"
        return values
