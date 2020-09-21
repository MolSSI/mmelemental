from .base import SimBase
from mmelemental.models.molecule.mm_molecule import Molecule
from pydantic import Field
from typing import Any, Tuple, Union

class Dynamics(Base):
    """ Molecular dynamics parameter input schema."""

    # Velocity fields
    gen_vel: Union[bool, Tuple[bool]]  = Field(None, description = 'Generate velocities from Maxwell distribution(s) at temperature(s) gen_temp, with random seed(s) gen_seed.')
    gen_temp: Union[float, Tuple[float]]  = Field(None, description = 'Temperature(s) in Kelvin for Maxwell distribution(s) used to generate gen_vel.')
    gen_seed: Union[int, Tuple[int]]  = Field(None, description = 'Used to initialize random generator for random velocities for each group.')


    # Thermostat fields
    temp: Tuple[float] = Field(None, description = 'Temperature(s) for thermostat(s) in Kelvin.')
    
    # Barostat fields
    press_couple: str = Field(None, description = 'Pressure coupling method.')
    press_couple_type: str = Field(None, description = 'Specifies the type of isotropy used for pressure coupling.')
    press_couple_freq: int = Field(None, description = 'Frequency used for coupling the pressure.')
    press_couple_relax: float = Field(None, description = 'The relaxation time constant in picoseconds used for pressure coupling (one value for all directions?).')
    compressibility: Union[float,Tuple[float]] = Field(None, description = 'Compressibility factor(s) in 1/bar. The number of required values is implied by press_couple_type.')
    press: Union[float,Tuple[float]] = Field(None, description = 'Pressure(s) used for barostate in bar. The number of required values is implied by press_couple_type.')
    # NAMD fields
    use_group_press: bool = Field(None, description = 'Pressure is calculated using either the atomic virial and kinetic energy (the default) or a hydrogen-group based pseudo-molecular virial and kinetic energy')
    use_flex_cell: bool = Field(None, description = 'Allows the three orthogonal dimensions of the periodic cell to fluctuate independently when this option is enabled.')
    use_const_ratio: bool = Field(None, description = 'When enabled, the ratio of unit cell in the xy plane is kept constant while allowing fluctuations along all axes.')
    use_const_area: bool = Field(None, description = 'When enabled, the ddimension of the unit cell in the xy is kept plane constant while allowing fluctuations along the z axis.')
    LangevinPistonTemp 
    SurfaceTensionTarget 
    StrainRate 
    ExcludeFromPressure 
    
    # GROMACS fields
    press_scale: str = Field(None, description = 'Which reference coordinates to scale the matrix of the pressure coupling with. Options: no, all, or com.')
