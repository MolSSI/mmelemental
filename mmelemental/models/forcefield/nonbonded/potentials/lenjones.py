from pydantic import Field, root_validator
from typing import Optional
import qcelemental
from mmelemental.models.base import ProtoModel

__all__ = ["LennardJones"]


class LennardJones(ProtoModel):

    epsilon: qcelemental.models.types.Array[float] = Field(
        ...,
        description="The epsilon (well depth) Lennard-Jones parameter. Default unit is kJ/mol.",
    )
    epsilon_units: Optional[str] = Field(
        "kJ/mol",
        description="Units for the Lennard-Jones epsilon (well depth) constant.",
    )
    sigma: qcelemental.models.types.Array[float] = Field(
        ...,
        description="The distance at which the Lennard-Jones potential is 0. Default unit is angstroms.",
    )
    sigma_units: Optional[str] = Field(
        "angstrom", description="Units for the Lennard-Jones sigma constant."
    )

    @root_validator(allow_reuse=True)
    def _valid_length(cls, values):
        assert len(values["epsilon"].shape) == 1, "epsilon must be a 1D array!"
        assert len(values["sigma"].shape) == 1, "sigma must be a 1D array!"
        assert len(values["epsilon"]) == len(
            values["sigma"]
        ), "epsilon and sigma must be of equal length!"
        return values

    def dict(self, *args, **kwargs):
        kwargs["exclude"] = {"provenance"}
        return super().dict(*args, **kwargs)
