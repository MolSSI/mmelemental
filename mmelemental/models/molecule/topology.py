import numpy
from typing import Optional
from pydantic import Field


# MM models
from mmelemental.models.base import ProtoModel, Provenance, provenance_stamp
from mmelemental.types import Array

from mmelemental.util.data import (
    NUMPY_INT,
    NUMPY_FLOAT,
    NUMPY_UNI,
)

__all__ = ["Topology"]


class Topology(ProtoModel):
    symbols: Optional[Array[NUMPY_UNI]] = Field(  # type: ignore
        None,
        description="An ordered (natom,) array-like object of particle symbols. The index of "
        "this attribute sets the order for all other per-particle setting like ``geometry`` and the first "
        "dimension of ``geometry``.",
    )
    atom_labels: Optional[Array[NUMPY_UNI]] = Field(  # type: ignore
        None,
        description="Additional per-atom labels as an array of strings. Typical use is in "
        "model conversions, such as Elemental <-> Molpro and not typically something which should be user "
        "assigned. See the ``comments`` field for general human-consumable text to affix to the molecule.",
    )
    atomic_numbers: Optional[Array[numpy.dtype(NUMPY_INT)]] = Field(  # type: ignore
        None,
        description="An optional ordered 1-D array-like object of atomic numbers of shape (nat,). Index "
        "matches the 0-indexed indices of all other per-atom settings like ``symbols`` and ``real``. "
        "Values are inferred from the ``symbols`` list if not explicitly set. "
        "Ghostedness should be indicated through ``real`` field, not zeros here.",
    )
    mass_numbers: Optional[Array[numpy.dtype(NUMPY_INT)]] = Field(  # type: ignore
        None,
        description="An optional ordered 1-D array-like object of atomic *mass* numbers of shape (nat). Index "
        "matches the 0-indexed indices of all other per-atom settings like ``symbols`` and ``real``. "
        "Values are inferred from the most common isotopes of the ``symbols`` list if not explicitly set. "
        "If single isotope not (yet) known for an atom, -1 is placeholder.",
    )
    masses: Optional[Array[numpy.dtype(NUMPY_FLOAT)]] = Field(  # type: ignore
        None,
        description="The ordered array of particle masses. Index order "
        "matches the 0-indexed indices of all other per-atom fields like ``symbols`` and ``real``. If "
        "this is not provided, the mass of each atom is inferred from its most common isotope. If this "
        "is provided, it must be the same length as ``symbols`` but can accept ``None`` entries for "
        "standard masses to infer from the same index in the ``symbols`` field.",
    )
    masses_units: Optional[str] = Field(  # type: ignore
        "amu",
        description="Units for atomic masses. Defaults to unified atomic mass unit.",
    )
    molecular_charge: Optional[float] = Field(  # type: ignore
        0.0,
        description="The net electrostatic charge of the molecule. Default unit is elementary charge.",
    )
    molecular_charge_units: Optional[str] = Field(  # type: ignore
        "e", description="Units for molecular charge. Defaults to elementary charge."
    )
    connectivity: Optional[Array[numpy.dtype(f"{NUMPY_INT}, {NUMPY_INT}, {NUMPY_FLOAT}")]] = Field(  # type: ignore
        None,
        description="A list of bonds within the molecule. Each entry is a tuple "
        "of ``(atom_index_A, atom_index_B, bond_order)`` where the ``atom_index`` "
        "matches the 0-indexed indices of all other per-atom settings like ``symbols`` and ``real``. "
        "Bonds may be freely reordered and inverted.",
    )
