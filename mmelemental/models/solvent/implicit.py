from mmelemental.models.base import Base
from pydantic import Field


_all__ = ["Solvent"]


class Solvent(Base):
    implicit: bool = Field(
        ...,
        description="Sets the solvent to be implicitly represented in a simulation.",
    )
