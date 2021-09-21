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

__all__ = ["EAM"]


class EAM(ProtoModel):

    embed: Array[float] = Field(
        ...,
        description="Array of embedding (discretized) energy terms. Default unit is kJ/mol.",
    )
    embed_units: Optional[str] = Field(
        "kJ/mol",
        description="Units for the embedding energy term.",
        dimensionality=MASS_DIM * LENGTH_DIM ** 2 / (TIME_DIM ** 2 * SUBS_DIM),
    )
    potential: Array[float] = Field(
        ...,
        description="Array of (discretized) pair potential interaction functions. Default unit is kJ/mol.",
    )
    pair_units: Optional[str] = Field(
        "kJ/mol",
        description="Units for the pair potential interaction term.",
        dimensionality=MASS_DIM * LENGTH_DIM ** 2 / (TIME_DIM ** 2 * SUBS_DIM),
    )
    density: Array[float] = Field(
        ..., description="Array of (discretized) atomic electron densities."
    )

    @root_validator(allow_reuse=True)
    def _valid_length(cls, values):
        assert len(values["embed"].shape) == 1, "embed must be a 1D array!"
        assert len(values["density"].shape) == 1, "density must be a 1D array!"
        assert len(values["pair"].shape) == 2, "pair must be a 2D array!"
        assert len(values["embed"]) == len(
            values["density"]
        ), "embed and density must be of equal length!"
        return values
