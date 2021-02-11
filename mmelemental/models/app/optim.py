from mmelemental.models.app.base import SimInput
from pydantic import Field
from typing import Tuple, Union


class OptimInput(SimInput):
    tol: float = Field(
        None,
        description="Tolerance used to indicate when the optimization scheme has converd.",
    )
