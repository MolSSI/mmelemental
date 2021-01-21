from mmelemental.models.base import Base
from pydantic import Field
from typing import Any


class ToolkitFF(Base):
    ff: Any = Field(..., description="Toolkit-specific force field object.")

    class Config(Base.Config):
        arbitrary_types_allowed = True

    @property
    def dtype(self):
        raise NotImplementedError
