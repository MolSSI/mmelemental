from pydantic import Field, root_validator
from typing import Optional
import qcelemental
from ..nb_params import NonBondedParams

__all__ = ["LennardJones"]


class LennardJones(NonBondedParams):

    epsilon: qcelemental.models.types.Array[float] = Field(
        ...,
        description="The epsilon (well depth) Lennard-Jones parameter. Default unit is kJ/mol.",
    )
    epsilon_units: Optional[str] = Field(
        "kJ/mol", description="The epsilon (well depth) Lennard-Jones unit."
    )
    sigma: qcelemental.models.types.Array[float] = Field(
        ...,
        description="The distance at which the Lennard-Jones potential is 0. Default unit is angstroms.",
    )
    sigma_units: Optional[str] = Field(
        "angstrom", description="The Lennard-Jones sigma unit."
    )

    @root_validator
    def _valid_length(cls, values):
        assert len(values["epsilon"].shape) == 1, "epsilon must be a 1D array!"
        assert len(values["sigma"].shape) == 1, "sigma must be a 1D array!"
        assert len(values["epsilon"]) == len(
            values["sigma"]
        ), "epsilon and sigma must be of equal length!"
        return values
