from mmelemental.models.base import Base
from pydantic import Field
from typing import Any

class ToolkitMolecule(Base):
    mol: Any = Field(..., description = 'Toolkit-specific molecule object.')
    
    class Config(Base.Config):
        arbitrary_types_allowed = True

    @property
    def dtype(self):
        raise NotImplementedError