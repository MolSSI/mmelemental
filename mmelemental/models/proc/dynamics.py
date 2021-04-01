from mmelemental.models.proc.base import ProcInput
from pydantic import Field
from qcelemental.models.types import Array

__all__ = ["DynamicsInput"]


class DynamicsInput(ProcInput):
    """ Molecular dynamics parameter input schema."""

    # Velocity fields
    gen_vel: Array[bool] = Field(
        None,
        description="Generate velocities from Maxwell distribution(s) at temperature(s) gen_temp, with random seed(s) gen_seed.",
    )
    gen_temp: Array[float] = Field(
        None,
        description="Temperature(s) in Kelvin for Maxwell distribution(s) used to generate gen_vel.",
    )
    gen_seed: Array[int] = Field(
        None,
        description="Used to initialize random generator for random velocities for each group.",
    )

    # Temperature coupling fields
    temp: Array[float] = Field(
        None, description="Temperature(s) for thermostat(s) in Kelvin."
    )
    temp_method: str

    # Pressure coupling fields
    press_method: str = Field(None, description="Pressure coupling method.")
    press: Array[float] = Field(
        None,
        description="Pressure(s) used for barostate in bar. The number of required values is implied by press_couple_type.",
    )
    press_couple_type: str = Field(
        None, description="Specifies the type of isotropy used for pressure coupling."
    )
    press_couple_freq: int = Field(
        None, description="Frequency used for coupling the pressure."
    )
    compressibility: Array[float] = Field(
        None,
        description="Compressibility factor(s) in 1/bar. The number of required values is implied by press_couple_type.",
    )

    # GROMACS fields
    press_scale: str = Field(
        None,
        description="Which reference coordinates to scale the matrix of the pressure coupling with. Options: no, all, or com.",
    )
