import numpy
from typing import Optional
from pydantic import Field, constr
import hashlib
import json


# MM models
from mmelemental.models.base import ProtoModel, Provenance
from cmselemental.types import Array

from mmelemental.util.data import (
    float_prep,
    NUMPY_INT,
    NUMPY_FLOAT,
    NUMPY_UNI,
    CHARGE_NOISE,
    MASS_NOISE,
)
from mmelemental.util.units import (
    MASS_DIM,
    TIME_DIM,
    CURRENT_DIM,
)

__all__ = ["Topology"]

mmschema_topology_default = "mmschema_topology"


class Topology(ProtoModel):
    """A topological representation of a Molecule in MM. In addition to connectivity, this model contains data for particle
    symbols/labels, masses, and net charge. Useful for creating :class:``Trajectory`` objects.
    """

    schema_name: constr(
        strip_whitespace=True, regex=mmschema_topology_default
    ) = Field(  # type: ignore
        mmschema_topology_default,
        description=(
            f"The MMSchema specification to which this model conforms. Explicitly fixed as {mmschema_topology_default}."
        ),
    )
    schema_version: int = Field(  # type: ignore
        1,
        description="The version number of ``schema_name`` to which this model conforms.",
    )
    name: Optional[str] = Field(  # type: ignore
        None,
        description="Common or human-readable name to assign to this molecule. This field can be arbitrary.",
    )
    comment: Optional[str] = Field(  # type: ignore
        None,
        description="Additional comments for this molecule. Intended for pure human/user consumption and clarity.",
    )
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
        dimensionality=MASS_DIM,
    )
    molecular_charge: Optional[float] = Field(  # type: ignore
        0.0,
        description="The net electrostatic charge of the molecule. Default unit is elementary charge.",
    )
    molecular_charge_units: Optional[str] = Field(  # type: ignore
        "e",
        description="Units for molecular charge. Defaults to elementary charge.",
        dimensionality=CURRENT_DIM * TIME_DIM,
    )
    connectivity: Optional[Array[numpy.dtype(f"{NUMPY_INT}, {NUMPY_INT}, {NUMPY_FLOAT}")]] = Field(  # type: ignore
        None,
        description="A list of bonds within the molecule. Each entry is a tuple "
        "of ``(atom_index_A, atom_index_B, bond_order)`` where the ``atom_index`` "
        "matches the 0-indexed indices of all other per-atom settings like ``symbols`` and ``real``. "
        "Bonds may be freely reordered and inverted.",
    )

    def __repr_args__(self) -> "ReprArgs":
        return [("name", self.name), ("hash", self.get_hash()[:7])]

    def __eq__(self, other):
        """
        Checks if two models are identical. This is a molecular identity defined
        by scientific terms, and not programing terms, so it's less rigorous than
        a programmatic equality or a memory equivalent `is`.
        """

        if isinstance(other, dict):
            other = self.__class__(**other)
        elif isinstance(other, self.__class__):
            pass
        else:
            raise TypeError(
                f"Comparison between {self.__class__} and {type(other)} is not supported."
            )

        return self.get_hash() == other.get_hash()

    @property
    def hash_fields(self):
        return [
            "symbols",
            "masses",
            "masses_units",
            "molecular_charge",
            "molecular_charge_units",
            "connectivity",
        ]

    def get_hash(self):
        """
        Returns the hash of the molecule.
        """

        m = hashlib.sha1()
        concat = ""

        # np.set_printoptions(precision=16)
        for field in self.hash_fields:
            data = getattr(self, field)
            if data is not None:
                if field == "molecular_charge":
                    data = float_prep(data, CHARGE_NOISE)
                elif field == "masses":
                    data = float_prep(data, MASS_NOISE)

                concat += json.dumps(data, default=lambda x: x.ravel().tolist())

        m.update(concat.encode("utf-8"))
        return m.hexdigest()
