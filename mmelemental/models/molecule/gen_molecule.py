from qcelemental import models
from pydantic import Field
from typing import Any

class ToolkitMolecule(models.ProtoModel):
    mol: Any = Field(..., description = 'Toolkit-specific molecule object.')
    dtype: str = Field(None, description = 'Data type for mol.')

    class Config(models.ProtoModel.Config):
        arbitrary_types_allowed = True

    @property
    def obj_type(self):
        return self.dtype