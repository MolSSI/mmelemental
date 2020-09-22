from mmelemental.models.base import Base
from mmelemental.models.sim.base import SimBase
from typing import List, Optional, Union, Dict
from pydantic import Field
from typing import Dict, List, Any, Optional, Tuple

class SimWriterInput(Base):
    model: SimBase = Field(..., description = 'Simulation schema model.')
    engine: Tuple[str, Any] = Field(..., description = 'Name of simulation engine and version: (engine_name, x.x.x).')
    filename: str = Field(..., description = 'Name of file to be written with the "to_file" method.')