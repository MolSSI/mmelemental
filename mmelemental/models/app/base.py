from mmelemental.models.base import Base
from mmelemental.models.molecule.mm_mol import Molecule
from mmelemental.models.solvent.implicit import Solvent
from mmelemental.models.forcefield import ForceField
from pydantic import Field
from typing import Tuple, Union, Dict


class SimInput(Base):
    """ Basic model for molecular simulation input parameters."""

    # System fields
    mol: Dict[str, Molecule] = Field(
        None,
        description="Molecular mechanics molecule object(s). See the :class:``Molecule`` class. "
        "Example: mol = {'ligand': Molecule, 'receptor': Molecule, 'solvent': Molecule}.",
    )
    cell: Tuple[Tuple[float], Tuple[float]] = Field(
        None,
        description="Cell dimensions in the form: ((xmin, ymin, ...), (xmax, ymax, ...))",
    )
    forcefield: Dict[str, ForceField] = Field(
        None, description='Forcefield object(s) for every Molecule defined in "mol".'
    )
    boundary: Tuple[str] = Field(
        None,
        description="Boundary conditions in all dimensions e.g. (periodic, periodic, periodic) imposes periodic boundaries in 3D.",
    )

    # Temporal fields
    march_method: str = Field(
        None,
        description="Time marching method to use e.g. velocity verlet, leapfrog, etc.",
    )
    nsteps: int = Field(None, description="Number of integration steps to perform.")
    step_size: float = Field(
        None, description="Step size e.g. timestep, minimization step size, etc."
    )

    # Bonds fields
    constraints: Union[Tuple[str], str] = Field(
        None, description="Groups of atoms/bonds/angles/etc. to constrain e.g. h-bonds."
    )
    constraints_method: Union[Tuple[str], str] = Field(
        None,
        description="Method used to constraint what's defined in 'constraints' e.g. LINCS or SHAKE.",
    )
    constraints_tol: float = Field(
        None, description="Tolerance used for constraint self-consistency."
    )

    # I/O fields
    traj: str = Field(
        None,
        description="Trajectory filename for coordinates and/or velocities, forces.",
    )
    traj_comp: str = Field(None, description="Compressed trajectory filename.")
    traj_vel: str = Field(
        None, description="Trajectory filename explicitly for velocities."
    )
    traj_for: str = Field(
        None, description="Trajectory filename explicitly for forces."
    )
    freq: int = Field(
        None, description="Frequency of writing to trajectory output file."
    )
    freq_vel: int = Field(
        None, description="Frequency of writing to velocity trajectory output file."
    )
    freq_for: int = Field(
        None, description="Frequency of writing to force trajectory output file."
    )
    energy: str = Field(None, description="Energy output filename.")
    freq_ene: int = Field(
        None, description="Frequency of writing to energy output file."
    )
    restart: str = Field(None, description="Restart output filename.")
    restart_freq: int = Field(
        None, description="Frequency of writing to restart output file."
    )
    restart_file: str = Field(None, description="Name of restart output file.")
