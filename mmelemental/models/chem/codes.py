from mmelemental.models.base import ProtoModel
from pydantic import Field, ValidationError
from mmelemental.util.decorators import require
from typing import Optional

__all__ = ["ChemCode"]


class ChemCode(ProtoModel):
    code: str = Field(
        ...,
        description="A chemical code that describes a molecule or molecular pattern e.g. smiles, smarts, etc. "
        "See the :class:``Identifiers`` class for supported codes.",
    )
    dtype: Optional[str] = Field(
        None, description="Data type e.g. smiles, smarts, etc. "
    )

    class _CodesSupported:
        codes = (
            "Smiles",
            "Smarts",
            "Inchi",
            "FASTA",
            "HELM",
            "Sequence",
        )  # must replace with Identifiers

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def __repr__(self):
        return self.code

    @require("rdkit")
    # @validator("code")
    def valid_code(cls, code):
        from rdkit import Chem

        for ccode in ChemCode._CodesSupported.codes:
            function = getattr(Chem, f"MolFrom{ccode}")
            if function(code):
                return code
        raise ValidationError

    @property
    @require("rdkit")
    def guess_dtype(self):
        from rdkit import Chem

        for ccode in ChemCode._CodesSupported.codes:
            function = getattr(Chem, f"MolFrom{ccode}")
            if function(self.code):
                return ccode
        raise ValidationError
