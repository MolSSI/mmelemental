from pydantic import Field, constr, validator
import importlib
import hashlib
import json
from typing import Any, List, Dict, Tuple, Optional, Union
import numpy
from pathlib import Path
from cmselemental.types import Array

from mmelemental.util.data import (
    float_prep,
    NUMPY_UNI,
    NUMPY_INT,
    NUMPY_FLOAT,
    GEOMETRY_NOISE,
    MASS_NOISE,
    CHARGE_NOISE,
)
from mmelemental.util.units import CURRENT_DIM, TIME_DIM, MASS_DIM

# MM models
from mmelemental.models.base import ProtoModel, Provenance, provenance_stamp
from mmelemental.models.util.output import FileOutput
from .nonbonded import NonBonded
from .bonded import Bonds, Angles, Dihedrals, DihedralsImproper

_trans_nfound_msg = "MMElemental translation requires mmic_translator. \
Solve by: pip install mmic_translator"

mmschema_forcefield_default = "mmschema_forcefield"


__all__ = ["ForceField"]


class ForceField(ProtoModel):
    schema_name: constr(
        strip_whitespace=True, regex="^(mmschema_forcefield)$"
    ) = Field(  # type: ignore
        mmschema_forcefield_default,
        description=(
            f"The MMSchema specification to which this model conforms. Explicitly fixed as {mmschema_forcefield_default}."
        ),
    )
    schema_version: int = Field(  # type: ignore
        1,
        description="The version number of ``schema_name`` to which this model conforms.",
    )
    author: Optional[str] = Field(  # type: ignore
        None,
        description="Author name to assign to the force field this model stores. This field can be arbitrary.",
    )
    name: Optional[str] = Field(  # type: ignore
        None,
        description="Common or human-readable name to assign to this model. This field can be arbitrary.",
    )
    version: Optional[str] = Field(  # type: ignore
        None,
        description="Version of the force field this model stores. This field can be arbitrary.",
    )
    comment: Optional[str] = Field(  # type: ignore
        None,
        description="Additional comments for this model. Intended for pure human/user consumption and clarity.",
    )
    symbols: Optional[Array[str]] = Field(  # type: ignore
        None,
        description="An ordered (natom,) list of particle (e.g. atomic elemental) symbols.",
    )
    nonbonded: Optional[Union[NonBonded, List[NonBonded]]] = Field(  # type: ignore
        None, description="Non-bonded parameters model."
    )
    bonds: Optional[Union[Bonds, List[Bonds]]] = Field(  # type: ignore
        None, description="2-body covalent bond model."
    )
    angles: Optional[Union[Angles, List[Angles]]] = Field(  # type: ignore
        None, description="3-body angular bond model."
    )
    dihedrals: Optional[Union[Dihedrals, List[Dihedrals]]] = Field(  # type: ignore
        None, description="4-body torsional bond model."
    )
    dihedrals_improper: Optional[Union[DihedralsImproper, List[DihedralsImproper]]] = Field(  # type: ignore
        None, description="Improper 4-body torsional bond model."
    )
    charges: Optional[Array[numpy.dtype(NUMPY_FLOAT)]] = Field(
        None, description="Atomic charges. Default unit is in elementary charge units."
    )
    charges_units: Optional[str] = Field(
        "e", description="Atomic charge unit.", dimensionality=CURRENT_DIM * TIME_DIM
    )
    masses: Optional[Array[numpy.dtype(NUMPY_FLOAT)]] = Field(  # type: ignore
        None,
        description="List of atomic masses. If not provided, the mass of each atom is inferred from its most common isotope. "
        "If this is provided, it must be the same length as ``symbols``.",
    )
    masses_units: Optional[str] = Field(  # type: ignore
        "amu",
        description="Units for atomic masses. Defaults to unified atomic mass unit.",
        dimensionality=MASS_DIM,
    )
    charge_groups: Optional[Array[numpy.dtype(NUMPY_INT)]] = Field(
        None, description="Charge groups per atom. Length of the array must be natoms."
    )
    exclusions: Optional[str] = Field(  # type: ignore
        None,
        description="Which pairs of bonded atoms to exclude from non-bonded calculations. \
    	The rules to apply in choosing bonded exclusions are specifed in the configuration file using the exclude parameter. The \
    	choices for exclusions are None, 1-2, 1-3, 1-4, etc. With None, no atom pairs are excluded. With 1-2, only atoms that are connected \
    	via a linear bond are excluded. With 1-3, any pair of atoms connected via a bond or bonded to a common third atom are excluded. \
    	With 1-4, any atoms that are connected by a linear bond, or a sequence of two bonds, or a sequence of three bonds, are excluded. \
    	With scaled1-4, exclusions are applied in the same manner as the 1-3 setting, but those pairs that are connected by a sequence of \
    	3 bonds are calculated using the modified 1-4 methods described rather than the standard force calculations.",
    )
    inclusions: Optional[str] = Field(  # type: ignore
        None,
        description="Which pairs of 1-4 excluded bonded atoms to include in non-bonded calculations.",
    )
    defs: Optional[List[str]] = Field(  # type: ignore
        None,
        description="Particle definition. For atomic forcefields, this could be the atom type (e.g. HH31) or SMIRKS (OFF) representation. "
        "The type names are associated with the atomic elements defined in other objects e.g. see the :class:``Molecule`` model.",
    )
    substructs: Optional[List[Tuple[str, int]]] = Field(
        None,
        description="A list of substructure names the particles belong to. E.g. [('ALA', 4), ('ACE', 0)] means atom1 belong to residue ALA (alanine) "
        "with residue number 4, while atom2 belongs to residue ACE (acetyl) with residue number 0.",
    )
    templates: Optional[Dict[str, List[str]]] = Field(
        None,
        description="A list of template definitions typically in terms of atom types. E.g. {'ACE': ['HH31', 'CH3', 'HH32', 'HH33', 'C', 'O']}.",
    )
    atomic_numbers: Optional[Array[numpy.int16]] = Field(  # type: ignore
        None,
        description="An optional ordered 1-D array-like object of atomic numbers of shape (nat,). Index "
        "matches the 0-indexed indices of all other per-atom settings like ``symbols``. "
        "Values are inferred from the ``symbols`` list if not explicitly set. ",
    )
    # Extras
    provenance: Provenance = Field(  # type: ignore
        provenance_stamp(__name__),
        description="The provenance information about how this object (and its attributes) were generated, "
        "provided, and manipulated.",
    )
    extras: Dict[str, Any] = Field(  # type: ignore
        None,
        description="Additional information to bundle with the object. Use for schema development and scratch space.",
    )

    class Config(ProtoModel.Config):
        repr_style = lambda self: [("name", self.name), ("hash", self.get_hash()[:7])]

    def __init__(self, **kwargs: Optional[Dict[str, Any]]) -> None:
        """
        Initializes the molecule object from dictionary-like values.
        Parameters
        ----------
        **kwargs : Any
            The values of the Molecule object attributes.
        """
        kwargs["schema_name"] = kwargs.pop("schema_name", "mmschema_forcefield")
        kwargs["schema_version"] = kwargs.pop("schema_version", 0)

        atomic_numbers = kwargs.get("atomic_numbers")
        if atomic_numbers is not None:
            if kwargs.get("symbols") is None:

                kwargs["symbols"] = [
                    qcelemental.periodictable.to_E(x) for x in atomic_numbers
                ]

        # We are pulling out the values *explicitly* so that the pydantic skip_defaults works as expected
        # All attributes set below are equivalent to the default set.
        super().__init__(**kwargs)

        values = self.__dict__

        if not values.get("name"):
            values["name"] = "forcefield"

    # Representation -> used by qcelemental's __repr__
    def __repr_args__(self) -> "ReprArgs":
        forms = [
            form.__class__.__name__
            for form in (
                self.nonbonded,
                self.bonds,
                self.angles,
                self.dihedrals,
                self.dihedrals_improper,
            )
            if form
        ]
        return [("name", self.name), ("form", forms), ("hash", self.get_hash()[:7])]

    # Validators
    @validator("charges")
    def _charges_length(cls, v, values):
        assert len(v.shape) == 1, "Atomic charges must be 1D array!"
        return v

    # Constructors
    @classmethod
    def from_file(
        cls,
        filename: str,
        dtype: Optional[str] = None,
        translator: Optional[str] = None,
        **kwargs,
    ) -> "ForceField":
        """
        Constructs a ForceField object from a file.

        Parameters
        ----------
        filename: str
            The topology or FF filename to build from.
        dtype: Optional[str], optional
            The type of file to interpret e.g. psf. If unset, mmelemental attempts to discover the file type.
        translator: Optional[str], optional
            Translator name e.g. mmic_parmed. Takes precedence over dtype. If unset, MMElemental attempts
            to find an appropriate translator if it is registered in the :class:`TransComponent` class.
        **kwargs: Optional[Dict[str, Any]], optional
            Any additional keywords to pass to the constructor.
        Returns
        -------
        ForceField
            A constructed ForceField object.

        """

        file_ext = Path(filename).suffix if filename else None

        if file_ext in [".json"]:

            dtype = file_ext.removeprefix(".")
            # Raw string type, read and pass through
            if dtype == "json":
                with open(filename, "r") as infile:
                    data = json.load(infile)
                dtype = "dict"
            else:
                raise KeyError(f"Data type not supported: {dtype}.")

            return cls.from_data(data, dtype=dtype, **kwargs)

        fileobj = FileOutput(path=filename) if filename else None

        dtype = dtype or fileobj.ext.removeprefix(".")
        ext = "." + dtype

        # Generic translator component
        try:
            from mmic_translator.components import TransComponent
        except Exception:
            TransComponent = None

        if not translator:
            if not TransComponent:
                raise ModuleNotFoundError(_trans_nfound_msg)
            from mmic_translator import reg_trans

            reg_trans = list(reg_trans)

            while not translator:
                translator = TransComponent.find_ffread_tk(ext, trans=reg_trans)
                if not translator:
                    raise ValueError(
                        f"Could not read top file with ext {ext}. Please install an appropriate translator."
                    )
                # Make sure we can import the translator module
                if importlib.util.find_spec(translator):
                    mod = importlib.import_module(translator)

        elif importlib.util.find_spec(translator):
            mod = importlib.import_module(translator)

        tkff_class = mod._classes_map.get("ForceField")

        if not tkff_class:
            raise ValueError(
                f"No ForceField model found while looking in translator: {translator}."
            )

        tkff = tkff_class.from_file(
            filename=fileobj.abs_path if fileobj else None,
            dtype=dtype,
        )

        return cls.from_data(tkff, dtype=tkff.dtype, **kwargs)

    @classmethod
    def from_data(cls, data: Any, **kwargs) -> "ForceField":
        """
        Constructs a ForceField object from a data object.

        Parameters
        ----------
        data: Any
            Data to construct ForceField from.
        **kwargs: Optional[Dict[str, Any]], optional
            Additional kwargs to pass to the constructors.

        Returns
        -------
        ForceField
            A constructed ForceField object.

        """
        if isinstance(data, dict):
            kwargs.pop("dtype", None)  # remove dtype if supplied
            kwargs.update(data)
            return cls(**kwargs)

        return data.to_schema(**kwargs)

    def to_file(
        self,
        filename: str,
        dtype: Optional[str] = None,
        translator: Optional[str] = None,
        **kwargs: Dict[str, Any],
    ) -> None:
        """Writes the ForceField to a file.

        Parameters
        ----------
        filename : str
            The filename to write to
        dtype : Optional[str], optional
            The type of file to write (e.g. psf, top, etc.), attempts to infer dtype from
            file extension if not provided.
        translator: Optional[str], optional
            Translator name e.g. mmic_parmed. Takes precedence over dtype. If unset, MMElemental attempts
            to find an appropriate translator if it is registered in the :class:`TransComponent` class.
        **kwargs: Optional[str, Dict], optional
            Additional kwargs to pass to the constructor.

        """
        if not dtype:
            from pathlib import Path

            ext = Path(filename).suffix
        else:
            ext = "." + dtype

        mode = kwargs.pop("mode", "w")

        if ext == ".json":
            stringified = self.json(**kwargs)

            with open(filename, mode) as fp:
                fp.write(stringified)

            return

        # Generic translator component
        try:
            from mmic_translator.components import TransComponent
        except Exception:
            TransComponent = None

        if not translator:
            if not TransComponent:
                raise ModuleNotFoundError(_trans_nfound_msg)
            translator = TransComponent.find_ffwrite_tk(ext)

        if not translator:
            raise NotImplementedError(
                f"File extension {ext} not supported with any installed translators."
            )

        tkff = self.to_data(translator=translator, **kwargs)
        tkff.to_file(filename, dtype=dtype, **kwargs)  # pass dtype?

    def to_data(
        self,
        dtype: Optional[str] = None,
        translator: Optional[str] = None,
        **kwargs: Dict[str, Any],
    ) -> "ToolkitModel":
        """
        Constructs a toolkit-specific forcefield from MMSchema ForceField.
        Which toolkit-specific component is called depends on which package is installed on the system.

        Parameters
        ----------
        translator: Optional[str], optional
            Translator name e.g. mmic_parmed. Takes precedence over dtype. If unset, MMElemental attempts
            to find an appropriate translator if it is registered in the :class:`TransComponent` class.
        dtype: str, optional
            Data type e.g. mdanalysis, parmed, etc.
        **kwargs: Optional[Dict[str, Any]]
            Additional kwargs to pass to the constructors.
        Results
        -------
        ToolkitModel
            Toolkit-specific ForceField object

        """
        try:
            from mmic_translator.components import TransComponent
        except Exception:
            TransComponent = None

        if not translator:
            if not dtype:
                raise ValueError(
                    f"Either translator or dtype must be supplied when calling {__name__}."
                )
            if not TransComponent:
                raise ModuleNotFoundError(_trans_nfound_msg)
            translator = TransComponent.find_trans(dtype)

        if importlib.util.find_spec(translator):
            mod = importlib.import_module(translator)
            tkff = mod._classes_map.get("ForceField")

            if not tkff:
                raise ValueError(
                    f"No ForceField model found while looking in translator: {translator}."
                )

            return tkff.from_schema(self)
        else:
            raise NotImplementedError(
                f"Translator {translator} not available. Make sure it is properly installed."
            )

    def __hash__(self):
        return hash(self.get_hash())

    def __eq__(self, other):
        """
        Checks if two models are identical. This is a forcefield identity defined
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
            "defs",
            "masses",
            "masses_units",
            "charges",
            "charges_units",
            "nonbonded",
            "bonds",
            "angles",
            "dihedrals",
            "dihedrals_improper",
            "exclusions",
            "inclusions",
        ]

    def get_hash(self):
        """
        Returns the hash of the force field object.
        """

        m = hashlib.sha1()
        concat = ""

        # np.set_printoptions(precision=16)
        for field in self.hash_fields:
            data = getattr(self, field)
            if data is not None:
                if field == "symbols" or field == "defs":
                    concat += json.dumps(data, default=lambda x: x.ravel().tolist())
                if field == "charges":
                    data = float_prep(data, CHARGE_NOISE)
                    concat += json.dumps(data, default=lambda x: x.ravel().tolist())
                if field == "masses":
                    data = float_prep(data, MASS_NOISE)
                    concat += json.dumps(data, default=lambda x: x.ravel().tolist())
                if field in (
                    "nonbonded",
                    "bonds",
                    "angles",
                    "dihedrals",
                    "dihedrals_improper",
                ):
                    if not isinstance(data, dict):
                        data = data.dict()
                concat += json.dumps(
                    data, default=lambda x: x.ravel().tolist()
                )  # if serialization fails, assume type is numpy.ndarray

        m.update(concat.encode("utf-8"))
        return m.hexdigest()

    @property
    def is_topology(self):
        """Returns True if model contains "topological" (i.e. subset of assigned ff params) data rather than forcefield definition."""
        return True if self.defs is None else False
