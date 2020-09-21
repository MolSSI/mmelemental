from mmelemental.models.base import Base
from mmelemental.models.molecule.mm_molecule import Molecule
from mmelemental.models.solvent.implicit import Solvent
from mmelemental.models.forcefield.atomtype import ForceField
from pydantic import Field
from typing import Any, Tuple, Union

class SimBase(Base):
    """ Basic model for molecular simulation input parameters."""

    # System fields
    mol: Union[Tuple[Molecule], Molecule] = Field(..., description = 'Molecular mechanics molecule objects.')
    cell: Tuple[float] = Field(None, description = 'Cell dimensions.')
    forcefield: Tuple[ForceField] = Field(None, description = 'Forcefield objects for every Molecule defined in "mol".')
    solvent: Tuple[Solvent] = Field(None, description = 'Solvent objects.')

    # Temporal fields
    march_method: str = Field(None, description = 'Time marching method to use e.g. velocity verlet, leapfrog, etc.') 
    nsteps: int = Field(None, description = "Number of integration steps to perform.")
    step_size: float = Field(None, description = 'Step size e.g. timestep, minimization step size, etc.')
    step_first: int = Field(None, description = 'The number of the first timestep typically used only when a simulation is a continuation of a previous simulation.')
    steps_cycle: int = Field(None, description = 'Number of timesteps in each cycle. Each cycle represents the number of timesteps between atom reassignments e.g. for NAMD.')

    # Bonds fields
    constraints: Union[Tuple[str], str] = Field(None, description = 'Groups of atoms/bonds/angles/etc. to constrain e.g. h-bonds.') 
    constraints_method: Union[Tuple[str], str] = Field(None, description = "Method used to constraint what's defined in 'constraints' e.g. LINCS or SHAKE.")
    constraints_tol: float = Field(None, description = "Tolerance used for constraint self-consistency.")
    constraints_order: int = Field(None, description = "Highest order in the expansion of the constraint coupling matrix.")
    constraints_angle: bool = Field(None, description = 'Maximum angle that a bond can rotate.')
    constraints_iter: int = Field(None, description = " Number of iterations to correct for rotational lengthening e.g. in LINCS.")
    constraints_start: bool = Field(None, description = 'If set to True, constraints are applied to the start configuration and reset shells. Set to False for continued simulations.')
    bond_potential: str = Field(None, description = 'Potential to be used for bonds e.g. morse, harmonic, etc.')

    # Algorithmic fields
    freq_list: int = Field(None, description = 'Frequency to update the neighbor list and long range forces.')
    cutoff_scheme: str = Field(None, description = 'Neighbor searching scheme e.g. Verlet.')
    ns_type: str = Field(None, description = 'Method to determine neighbor list e.g. grid.')
    longrange_type: str = Field(None, description = 'Method for computing long range electrostatic interactions e.g. PME.')
    cutoff_coulomb: float = Field(None, description = 'Short-range electrostatic cut-off.')
    cutoff_vwd: float = Field(None, description = 'Short-range Van der Waals cut-off.')
    boundary: Tuple[str] = Field(None, description = 'Boundary conditions in all dimensions e.g. (pbc, pbc, pbc) imposed periodic boundaries in 3D.')

    # I/O fields
    traj: str = Field(None, description = "Trajectory filename for coordinates and/or velocities, forces.")
    traj_comp: str = Field(None, description = "Compressed trajectory filename.")
    traj_vel: str = Field(None, description = "Trajectory filename explicitly for velocities.")
    traj_for: str = Field(None, description = "Trajectory filename explicitly for forces.")
    freq: int = Field(None, description = "Frequency of writing to trajectory output file.")
    freq_vel: int = Field(None, description = "Frequency of writing to velocity trajectory output file.")
    freq_for: int = Field(None, description = "Frequency of writing to force trajectory output file.")
    energy: str = Field(None, description = "Energy filename.")
    freq_ene: int = Field(None, description = "Frequency of writing to energy output file.")
    restart: str = Field(None, description = "Restart filename.")
    restart_freq: int = Field(None, description = "Frequency of writing to restart output file.")
    restart_file: str = Field(None, description = "Name of restart output file.")