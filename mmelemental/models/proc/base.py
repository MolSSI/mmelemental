from mmelemental.models.base import ProtoModel
from mmelemental.models.molecule.mm_mol import Molecule
from mmelemental.models.collect.sm_ensem import Ensemble
from mmelemental.models.collect.mm_traj import Trajectory
from mmelemental.models.forcefield import ForceField
from pydantic import Field, constr
from typing import List, Union, Dict, Optional, Any

__all__ = ["ProcInput", "ProcOutput"]

mmschema_proc_input_default = "mmschema_proc_input"


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
    schema_name: constr(
        strip_whitespace=True, regex="^(mmschema_proc_input)$"
    ) = Field(  # type: ignore
        mmschema_proc_input_default,
        description=(
            f"The MMSchema specification to which this model conforms. Explicitly fixed as {mmschema_proc_input_default}."
        ),
    )
    schema_version: int = Field(  # type: ignore
        0,
        description="The version number of ``schema_name`` to which this model conforms.",
    )
    kwargs: Optional[Dict[str, Any]] = Field(
        None, description="Additional keyword arguments to pass to the constructors."
    )


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
    warnings: Optional[List[str]] = Field(
        None, description="Warning messages generated from the conversion."
    )
    stdout: str = Field(..., description="Standard output.")
    stderr: Optional[str] = Field(None, description="Standard error.")
    log: Optional[str] = Field(None, description="Logging output.")
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
