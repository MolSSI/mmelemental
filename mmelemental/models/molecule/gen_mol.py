from mmelemental.models.base import Base
from pydantic import Field
from typing import Any
import importlib
import inspect


class ToolkitMol(Base):
    """ An abstract base class that acts as a wrapper for toolkit molecules """

    mol: Any = Field(
        ..., description="Toolkit-specific molecule object."
    )  # Will be validated during runtime

    class Config(Base.Config):
        arbitrary_types_allowed = True

    @property
    def toolkit(self):
        return type(self.mol).__module__

    @property
    def translator(self):
        name, _ = self.__module__.split(".", 1)
        return name

    @property
    def path(self):
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
        if isinstance(self.mol, self.dtype):
            return self.mol
        raise ValidationError

    @staticmethod
    def check_name(name) -> str:
        """ Returns atom name of langth 4 characters. """
        assert len(name) <= 4

        if len(name) != 4:
            if len(name) == 1:
                name = " " + name + "  "
            elif len(name) == 2:
                name = " " + name + " "
            elif len(name) == 3:
                name = " " + name
        return name
