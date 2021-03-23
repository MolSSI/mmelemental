from mmelemental.models.base import ProtoModel
from mmelemental.models.molecule.mm_mol import Molecule
from mmelemental.models.collect.sm_ensem import Ensemble
from mmelemental.models.collect.mm_traj import Trajectory, TrajInput
from mmelemental.models.solvent.implicit import Solvent
from mmelemental.models.forcefield import ForceField
from pydantic import Field, validator
from typing import Tuple, List, Union, Dict, Optional, Any

__all__ = ["ProcInput", "ProcOutput"]


class ProcInput(ProtoModel):
    """ Basic model for molecular simulation input parameters."""

    # Generic fields
    engine: Optional[str] = Field(
        None,
        description="Engine name to use in the procedure e.g. OpenMM.",
    )
    engine_version: Optional[str] = Field(
        None, description="Supported engine version. e.g. >= 3.4.0."
    )
    component: Optional[str] = Field(
        None,
        description="Component name to use in the procedure e.g. mmic_openmm.",
    )
    schema_version: Optional[str] = Field(
        None, description="Supported schema version. e.g. >= 1.2.0."
    )  # we need this? yah, but non-MMSchemas
    schema_name: Optional[str] = Field("MMSchema", description="Schema name.")
    kwargs: Optional[Dict[str, Any]] = Field(
        None, description="Additional keyword arguments to pass to the constructors."
    )

    # System fields
    mol: Dict[str, Molecule] = Field(
        None,
        description="Molecular mechanics molecule object(s). See the :class:``Molecule`` class. "
        "Example: mol = {'ligand': Molecule, 'receptor': Molecule, 'solvent': Molecule}.",
    )
    forcefield: Optional[Union[Dict[str, ForceField], Dict[str, str]]] = Field(
        None,
        description='Forcefield object(s) or name(s) for every Molecule defined in "mol".',
    )
    cell: Optional[Tuple[Tuple[float], Tuple[float]]] = Field(
        None,
        description="Cell dimensions in the form: ((xmin, ymin, ...), (xmax, ymax, ...))",
    )
    boundary: Tuple[str] = Field(
        None,
        description="Boundary conditions in all dimensions e.g. (periodic, periodic, periodic) imposes periodic boundaries in 3D.",
    )

    # I/O fields
    trajectory: Optional[Dict[str, TrajInput]] = Field(
        None,
        description="Trajectories to write for quantity 'key' every 'value' steps. E.g. {'geometry': 10, 'velocities': 100, 'forces': 50} "
        "produces 3 trajectory objects storing positions every 10 steps, velocities, every 100 steps, and forces every 50 steps. A comma "
        "seperator is used to indicate multiple variables stored in the same trajectory e.g. {'geometry,y,z': 1} produces a single trajectory object which stores the x, y, and z positions "
        "every step.",
    )

    @validator("forcefield")
    def _valid_ff(cls, v, values, **kwargs):
        for name in values["mol"]:
            if name not in v:
                raise ValueError(f"{name} does not have a defined force field.")
        assert len(v) == len(values["mol"]), (
            "Every molecule should have a single force field definition. "
            + f"{len(values['mol'])} molecules defined using {len(v)} force fields."
        )

        return v


class ProcOutput(ProtoModel):
    """ Basic model for molecular simulation output."""

    component: str = Field(
        None,
        description="Component name used in the procedure e.g. mmic_openmm.",
    )
    engine: Optional[str] = Field(
        None,
        description="Engine name used in the procedure e.g. OpenMM.",
    )
    engine_version: Optional[str] = Field(
        None, description="Engine version used in the procedure e.g. >= 3.4.0."
    )
    warning: Optional[List[str]] = Field(
        None, description="Warning messages generated from the conversion."
    )

    mol: Optional[Dict[str, Molecule]] = Field(
        None,
        description="Molecular mechanics molecule object(s). See the :class:``Molecule`` class. "
        "Example: mol = {'ligand': Molecule, 'receptor': Molecule, 'solvent': Molecule}.",
    )
    forcefield: Optional[Union[Dict[str, ForceField], Dict[str, str]]] = Field(
        None,
        description='Forcefield object(s) or name(s) for every Molecule defined in "mol".',
    )
    ensemble: Optional[Dict[str, Ensemble]] = Field(
        None,
        description="Ensemble output for a series of microstates of molecules. "
        "See the :class:``Ensemble`` class.",
    )
    trajectory: Optional[Dict[str, Trajectory]] = Field(
        None,
        description="Trajectory output representing a series of snapshots of the system at "
        "different timesteps. See the :class:``Trajectory`` class.",
    )
    trajectory_units: Optional[Dict[str, str]] = Field(
        None,
        description="Trajectory units. Any unit supported by pint is allowed.",
    )
    observable: Optional[Dict[str, List[float]]] = Field(
        None,
        description="Stores any observable or physical variable not accounted for in the schema. "
        "e.g. ligand scores used in docking simulations.",
    )
    observable_units: Optional[Dict[str, str]] = Field(
        None,
        description="Observable units. Any unit supported by pint is allowed.",
    )
