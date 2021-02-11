from qcelemental import models
from pydantic import Field, ValidationError
from typing import Dict, Optional
from mmelemental.extras import get_information
from typing import Optional, Any
import importlib
import inspect


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


class Base(models.ProtoModel):
    provenance: Provenance = Field(
        provenance_stamp(__name__),
        description="The provenance information about how this object (and its attributes) were generated, "
        "provided, and manipulated.",
    )

    class Config(models.ProtoModel.Config):
        canonical_repr = True
        extra = "allow"


class Nothing(Base):
    ...


class ToolkitModel(Base):
    """ An abstract base class that acts as a wrapper for toolkit molecules """

    # class Config(Base.Config):
    #    arbitrary_types_allowed = True

    data: Any = Field(
        ..., description="Toolkit-specific trajectory object."
    )  # validator added in subclasses
    units: Optional[Dict] = Field(
        None, description="Unit system for the stored physical properties in data."
    )

    @property
    def toolkit(self) -> str:
        """ Returns the path module that defines the data type object. """
        return type(self.data).__module__

    @property
    def translator(self) -> str:
        name, _ = self.__module__.split(".", 1)
        return name

    @property
    def path(self) -> str:
        return self.__module__ + "." + self.__name__

    @property
    def components(self):
        comp_mod = importlib.import_module(self.translator + ".components")
        return inspect.getmembers(comp_mod, inspect.isclass)

    @property
    def models(self):
        mod = importlib.import_module(self.translator + ".models")
        return inspect.getmembers(mod, inspect.isclass)

    def check_type(self):
        if isinstance(self.data, self.dtype):
            return self.data
        raise ValidationError
