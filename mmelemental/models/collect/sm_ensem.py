from pydantic import Field
from typing import Union, Optional, List, Dict, Any
from qcelemental.models.types import Array
from mmelemental.models.molecule.mm_mol import Molecule
from mmelemental.models.base import ProtoModel


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
    pot_energy: Optional[Array[float]] = Field(
        None,
        description="Total system potential energy. Default unit is KiloJoules/mol.",
    )
    pot_energy_units: Optional[str] = Field(
        "kJ/mol", description="Potential energy units. Defaults to KiloJoules/mol."
    )
    observables: Optional[Dict[str, Any]] = Field(
        None,
        description="Observables or physical variables not accounted in the schema \
            e.g. ligand scores used in docking simulations.",
    )


class Ensemble(ProtoModel):
    mol: Optional[Union[List[Molecule], Molecule]] = Field(
        None,
        description="Single or multiple :class:``Molecule`` object(s) representing the molecular topology.",
    )
    states: Optional[Dict[str, List[Microstate]]] = Field(
        None,
        description="Similar to Molecule but without the connectivity. Provides improved efficiency over the \
            latter. See :class:``Microstate``.",
    )
