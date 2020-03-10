from qcelemental import models
from pydantic import Field
from typing import Any

class ToolkitMolecule(models.ProtoModel):
    mol: Any = Field(None, description = 'toolkit-specific molecule object.')

    class Config(models.ProtoModel.Config):
        arbitrary_types_allowed = True