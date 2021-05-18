from pydantic import Field, root_validator
from typing import Optional
import qcelemental
from mmelemental.models.base import ProtoModel

__all__ = ["EAM"]


class EAM(ProtoModel):

    embed: qcelemental.models.types.Array[float] = Field(
        ...,
        description="Array of embedding (discretized) energy terms. Default unit is kJ/mol.",
    )
    embed_units: Optional[str] = Field(
        "kJ/mol", description="Units for the embedding energy term."
    )
    potential: qcelemental.models.types.Array[float] = Field(
        ...,
        description="Array of (discretized) pair potential interaction functions. Default unit is kJ/mol.",
    )
    pair_units: Optional[str] = Field(
        "kJ/mol", description="Units for the pair potential interaction term."
    )
    density: qcelemental.models.types.Array[float] = Field(
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

    def dict(self, *args, **kwargs):
        kwargs["exclude"] = {"provenance"}
        return super().dict(*args, **kwargs)
