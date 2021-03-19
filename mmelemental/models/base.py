from qcelemental import models
from pydantic import Field, ValidationError, validator
from typing import Dict, Optional
from mmelemental.extras import get_information
from typing import Optional


class Provenance(models.ProtoModel):
    """
    Provenance information.
    """

    creator: str = Field(..., description="The creator of the object.")
    version: Optional[str] = Field(None, description="The version of the creator.")
    routine: Optional[str] = Field(None, description="The routine of the creator.")

    class Config(models.ProtoModel.Config):
        canonical_repr = True
        extra = "allow"


def provenance_stamp(routine: str) -> Dict[str, str]:
    """Return dictionary satisfying QCSchema,
    https://github.com/MolSSI/QCSchema/blob/master/qcschema/dev/definitions.py#L23-L41
    with QCElemental's credentials for creator and version. The
    generating routine's name is passed in through `routine`.
    """
    return {
        "creator": "MMElemental",
        "version": get_information("version"),
        "routine": routine,
    }


class ProtoModel(models.ProtoModel):
    provenance: Provenance = Field(
        provenance_stamp(__name__),
        description="The provenance information about how this object (and its attributes) were generated, "
        "provided, and manipulated.",
    )

    def dict(self, *args, **kwargs):
        kwargs["by_alias"] = True
        kwargs["exclude_unset"] = False
        kwargs["exclude_none"] = True
        return super().dict(*args, **kwargs)

    @classmethod
    def get_units(cls):
        """ Returns model default units i.e. any Field name ending with _units. """
        return {
            val.name: val.default
            for key, val in cls.__fields__.items()
            if val.name.endswith("_units")
        }
