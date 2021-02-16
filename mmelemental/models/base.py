from qcelemental import models
from pydantic import Field, ValidationError, validator
from typing import Dict, Optional
from mmelemental.extras import get_information
from typing import Optional, Any
import importlib
import inspect
import abc


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
        # extra = "allow"


class Nothing(Base):
    ...


class ToolkitModel(Base, abc.ABC):
    """ An abstract base class that acts as a wrapper for toolkit molecules """

    data: Any = Field(
        ..., description="Toolkit-specific trajectory object."
    )  # validator added in subclasses
    units: Optional[Dict] = Field(
        None, description="Unit system for the stored physical properties in data."
    )

    @property
    @abc.abstractmethod
    def dtype(self):
        """ Returns the fundamental molecule object type. """
        raise NotImplementedError

    @classmethod
    @abc.abstractmethod
    def from_file(cls, filename: str = None, dtype: str = None, **kwargs):
        """ Constructs a data object from file(s). """
        raise NotImplementedError

    @classmethod
    @abc.abstractmethod
    def from_schema(cls, data: Any, version: Optional[str] = None, **kwargs):
        """ Constructs data object from MMSchema. """
        raise NotImplementedError

    @abc.abstractmethod
    def to_file(self, filename: str, dtype: str = None, **kwargs):
        """Writes the data object to a file.
        Parameters
        ----------
        filename : str
            The filename to write to
        dtype : Optional[str], optional
            File format
        **kwargs
            Additional kwargs to pass to the constructors.
        """
        raise NotImplementedError

    @abc.abstractmethod
    def to_schema(self, version: Optional[str] = None, **kwargs):
        """Converts the data object to MMSchema compliant object.
        Parameters
        ----------
        version: str, optional
            Schema specification version to comply with e.g. 1.0.1.
        **kwargs
            Additional kwargs to pass to the constructor.
        """
        raise NotImplementedError

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

    # @validator("data")
    def check_type(cls, data):
        if isinstance(data, cls.dtype):
            return data
        raise ValidationError
