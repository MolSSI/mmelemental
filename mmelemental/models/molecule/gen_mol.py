from mmelemental.models.base import Base
from pydantic import Field
from typing import Any

class ToolkitMol(Base):
    mol: Any = Field(..., description = 'Toolkit-specific molecule object.')
    
    class Config(Base.Config):
        arbitrary_types_allowed = True

    @property
    def dtype(self):
        raise NotImplementedError

    @staticmethod
    def check_name(name) -> str:
        """ Returns atom name of langth 4 characters. """
        assert len(name) <= 4

        if len(name) != 4:
            if len(name) == 1:
                name = ' ' + name + '  '
            elif len(name) == 2:
                name = ' ' + name + ' '
            elif len(name) == 3:
                name = ' ' + name
        return name