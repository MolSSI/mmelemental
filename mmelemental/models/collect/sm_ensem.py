from pydantic import Field
from typing import Optional, List, Dict
from qcelemental.models.types import Array
from mmelemental.models.molecule.mm_mol import Molecule
from mmelemental.models.base import ProtoModel, Provenance, provenance_stamp
import functools

__all__ = ["Microstate", "Ensemble"]


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
    mol: Optional[Dict[str, List[Molecule]]] = Field(
        None,
        description="Single or multiple :class:``Molecule`` object(s) representing the molecular topology.",
    )
    states: Optional[Dict[str, List[Microstate]]] = Field(
        None,
        description="Similar to Molecule but without the connectivity. Provides improved efficiency over the \
            latter. See :class:``Microstate``.",
    )
    provenance: Provenance = Field(
        default_factory=functools.partial(provenance_stamp, __name__),
        description="The provenance information about how this object (and its attributes) were generated, "
        "provided, and manipulated.",
    )
