from pydantic import Field, constr
from typing import Optional, List, Dict, Any
from qcelemental.models.types import Array
from mmelemental.models.molecule.mm_mol import Molecule
from mmelemental.models.base import ProtoModel, Provenance, provenance_stamp
import functools

__all__ = ["Microstate", "Ensemble"]

mmschema_ensemble_default = "mmschema_ensemble"


class Microstate(ProtoModel):
    geometry: Optional[Array[float]] = Field(
        ..., description="Atomic positions of length natoms. Default unit is Angstroms."
    )
    geometry_units: Optional[str] = Field(
        "angstrom", description="Units for atomic geometry. Defaults to Angstroms."
    )
    velocities: Optional[Array[float]] = Field(
        None,
        description="Atomic velocities of length natoms. Default unit is Angstroms/femotoseconds.",
    )
    velocities_units: Optional[str] = Field(
        "angstrom/fs",
        description="Units for atomic velocities. Defaults to Angstroms/femtoseconds.",
    )
    forces: Optional[Array[float]] = Field(
        None, description="Atomic forces of length natoms. KiloJoules/mol.Angstroms."
    )
    forces_units: Optional[str] = Field(
        "kJ/(mol*angstrom)",
        description="Units for atomic forces. Defaults to KiloJoules/mol.Angstroms",
    )


class Ensemble(ProtoModel):
    """A representation of ensembles in statistical mechanics. Useful for storing output
    from molecular docking, coarse-graining, etc.
    """

    schema_name: constr(
        strip_whitespace=True, regex=mmschema_ensemble_default
    ) = Field(  # type: ignore
        mmschema_ensemble_default,
        description=(
            f"The MMSchema specification to which this model conforms. Explicitly fixed as {mmschema_ensemble_default}."
        ),
    )
    schema_version: int = Field(  # type: ignore
        0,
        description="The version number of ``schema_name`` to which this model conforms.",
    )
    mols: Optional[Dict[str, List[Molecule]]] = Field(
        None,
        description="Single or multiple :class:``Molecule`` object(s) representing the molecular topology.",
        # + Molecule.__doc__,
    )
    states: Optional[Dict[str, List[Microstate]]] = Field(
        None,
        description="Similar to Molecule but without the connectivity. Provides improved efficiency over the \
            latter. See :class:``Microstate``.",
    )
    scores: List[float] = Field(
        ...,
        description="A list of scores for each state. Length must be equal to the number of states or mols.",
    )
    scores_units: Optional[str] = Field(None, description="Score function unit.")
    provenance: Provenance = Field(
        provenance_stamp(__name__),
        description="The provenance information about how this object (and its attributes) were generated, "
        "provided, and manipulated.",
    )
    extras: Dict[str, Any] = Field(  # type: ignore
        None,
        description="Additional information to bundle with this object. Use for schema development and scratch space.",
    )
