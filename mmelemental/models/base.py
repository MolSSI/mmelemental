from cmselemental import models
from pydantic import Field, root_validator
from typing import Dict, Optional
from mmelemental.extras import get_information
from cmselemental.util.decorators import classproperty
from typing import Optional, Any, Dict
import pint

__all__ = ["Provenance", "ProtoModel"]


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
    class Config:
        json_encoders: Dict[str, Any] = {pint.unit.UnitsContainer: lambda v: dict(v)}

    def dict(self, *args, **kwargs):
        kwargs["by_alias"] = True
        kwargs["exclude_unset"] = False
        kwargs["exclude_none"] = True
        return super().dict(*args, **kwargs)

    @classproperty
    def default_units(cls):
        """Returns class default units i.e. any Field name ending with _units."""
        return {
            val.name: val.default
            for val in cls.__fields__.values()
            if val.name.endswith("_units")
        }

    @property
    def units(self):
        """Returns instance (object) units i.e. any Field name ending with _units."""
        return {
            key: val for key, val in self.__dict__.items() if key.endswith("_units")
        }

    @root_validator
    def _valid_unit(cls, values: Dict[str, Any]) -> Dict[str, Any]:
        """Ensures all model units are physically valid and compatible with Pint. Returns
        field values with pint-compatible unit names."""
        schema_props = cls.schema()["properties"]
        for key, val in schema_props.items():
            if key.endswith("_units"):
                units = values.get(key)
                if units is not None:  # Make sure unit name/symbol has been supplied
                    quant = pint.Quantity(
                        1.0, units=units
                    )  # ensure pint supports this unit
                    dim = val.get("dimensionality")  # dim metadata must be defined
                    if dim is None:
                        raise AttributeError(
                            f"{key} Field does not store any dimensionality metadata."
                        )
                    unit_str = str(quant.u)
                    assert quant.check(
                        dim
                    ), f"Unit supplied ({unit_str}) for {key} has invalid dimensions {dim}."
                    values[key] = unit_str
        return values
