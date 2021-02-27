from pydantic import Field, root_validator
from typing import Optional
import qcelemental
from ..nb_params import NonBondedParams

__all__ = ["EAM"]


class EAM(NonBondedParams):

    embed: qcelemental.models.types.Array[float] = Field(
        ..., description="Embedding energy term. Default unit is kJ/mol."
    )
    embed_units: Optional[str] = Field(
        "kJ/mol", description="Units for the embedding energy term."
    )
    pair: qcelemental.models.types.Array[float] = Field(
        ..., description="Pair potential interaction. Default unit is kJ/mol."
    )
    pair_units: Optional[str] = Field(
        "kJ/mol", description="Units for the pair potential interaction."
    )
    density: qcelemental.models.types.Array[float] = Field(
        ..., description="Atomic electron density."
    )

    @root_validator
    def _valid_length(cls, values):
        assert len(values["embed"].shape) == 1, "embed must be a 1D array!"
        assert len(values["density"].shape) == 1, "density must be a 1D array!"
        assert len(values["pair"].shape) == 2, "pair must be a 2D array!"
        assert len(values["embed"]) == len(
            values["density"]
        ), "embed and density must be of equal length!"
        return values
