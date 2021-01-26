from mmelemental.models.base import Base
from mmelemental.models.input.base import SimInput
from pydantic import Field
from typing import Any, Tuple
import os


class SimWriterInput(Base):
    model: SimInput = Field(..., description="Simulation schema model.")
    engine: Tuple[str, Any] = Field(
        ..., description="Name of simulation engine and version: (engine_name, x.x.x)."
    )
    filename: str = Field(
        ...,
        description='Name of schema file to be written via "to_file" method. The base name is \
    						also used for coordinates, topology, etc. files.',
    )

    @property
    def filepath(self):
        return os.path.abspath(self.filename)
