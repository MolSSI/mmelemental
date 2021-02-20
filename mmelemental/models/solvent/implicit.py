from mmelemental.models.base import ProtoModel
from pydantic import Field


__all__ = ["Solvent"]


class Solvent(ProtoModel):
    implicit: bool = Field(
        ...,
        description="Sets the solvent to be implicitly represented in a simulation.",
    )
