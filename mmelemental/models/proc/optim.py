from mmelemental.models.proc.base import ProcInput
from mmelemental.models.base import ProtoModel
from pydantic import Field
from typing import Tuple, Union

class Frame(ProtoModel):
    pot_energy: Optional[List[float]] = Field(
        None,
        description="Total system potential energy. Default unit is KiloJoules/mol.",
    )


class FrameUnits(ProtoModel):
    pot_energy_units: Optional[str] = Field(
        "kJ/mol", description="Potential energy units. Defaults to KiloJoules/mol."
    )

class OptimInput(ProcInput):
    tol: float = Field(
        None,
        description="Tolerance used to indicate when the optimization scheme has converged.",
    )
    tol_units: float = Field(
        None,
        description="Tolerance unit e.g. kJ/mol."
    )
    observables: Optional[Frame] = Field(
        None,
        description="Observables or physical variables not accounted for in the schema. "
        "e.g. ligand scores used in docking simulations.",
    )
    observables_units: Optional[FrameUnits] = Field(
        None,
        description="Units observables. Any unit supported by pint is allowed.",
    )