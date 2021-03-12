from mmelemental.models.base import ProtoModel
from mmelemental.models.molecule.mm_mol import Molecule
from mmelemental.models.collect.sm_ensem import Ensemble
from mmelemental.models.collect.mm_traj import Trajectory
from mmelemental.models.solvent.implicit import Solvent
from mmelemental.models.forcefield import ForceField
from pydantic import Field
from typing import Tuple, List, Union, Dict, Optional

__all__ = ["SimInput", "SimOutput"]


class SimInput(ProtoModel):
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

    # Constraint fields
    bond_const: Optional[Dict[str, List[int]]] = Field(
        None,
        description="Specifies which bonds/angles/etc. in a molecule are constrained specified by their indices. E.g bond_const = {'solvent': [0,2,6]}.",
    )
    bond_const_method: Optional[str] = Field(
        None,
        description="Method used to constraint what's defined in 'constraints' e.g. LINCS or SHAKE.",
    )
    bond_const_tol: Optional[float] = Field(
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


class SimOutput(ProtoModel):
    """ Basic model for molecular simulation output."""

    simInput: SimInput = Field(
        ..., description="Simulation input used to generate the output."
    )
    ensemble: Ensemble = Field(
        None,
        description="Ensemble output for a series of microstates of molecules. "
        "See the :class:``Ensemble`` class.",
    )
    trajectory: Trajectory = Field(
        None,
        description="Trajectory output representing a series of snapshots of the system at "
        "different timesteps. See the :class:``Trajectory`` class.",
    )
    observables: Optional[Dict[str, List[float]]] = Field(
        None,
        description="Physical observables such as RMSD, energy, etc. E.g. observables={'RMSD':[...]}.",
    )
    observables_units: Optional[Dict[str, str]] = Field(
        None,
        description="Physical observables units. E.g. observables_units={'RMSD':'angstrom'}.",
    )
    pot_energy: Optional[List[float]] = Field(
        None,
        description="Total system potential energy. Default unit is KiloJoules/mol.",
    )
    pot_energy_units: Optional[str] = Field(
        "kJ/mol", description="Potential energy units. Defaults to KiloJoules/mol."
    )
    observables: Optional[Dict[str, List[float]]] = Field(
        None,
        description="Observables or physical variables not accounted for in the schema. "
        "e.g. ligand scores used in docking simulations.",
    )
    observables_units: Optional[Dict[str, str]] = Field(
        None,
        description="Units observables. Any unit supported by pint is allowed.",
    )
