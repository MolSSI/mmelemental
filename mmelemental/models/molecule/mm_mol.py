import qcelemental
import numpy
from typing import List, Tuple, Optional, Any, Dict, Union
from pydantic import Field, constr, validator
import importlib
from pathlib import Path
import hashlib
import json
import functools

# MM models
from mmelemental.models.util.output import FileOutput
from mmelemental.models.chem.codes import ChemCode
from mmelemental.models.base import Provenance, provenance_stamp, ProtoModel
from mmelemental.types import Array

__all__ = ["Molecule"]

# Rounding quantities for hashing
GEOMETRY_NOISE = 8
VELOCITY_NOISE = 8
MASS_NOISE = 6
CHARGE_NOISE = 4

_trans_nfound_msg = "MMElemental translation requires mmic & mmic_translator. \
Solve by: pip install mmic mmic_translator"
mmschema_molecule_default = "mmschema_molecule"


class Identifiers(qcelemental.models.molecule.Identifiers):
    """
    An extension of the qcelemental.models.molecule.Identifiers for RDKit constructors.
    See `link <https://rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html>`_ for more info.
    """

    smiles: Optional[Union[ChemCode, str]] = Field(
        None, description="A simplified molecular-input line-entry system code."
    )
    smarts: Optional[Union[ChemCode, str]] = Field(
        None,
        description="A SMILES arbitrary target specification code for defining substructures.",
    )
    inchi: Optional[Union[ChemCode, str]] = Field(
        None, description="An international chemical identifier code."
    )
    sequence: Optional[Union[ChemCode, str]] = Field(
        None,
        description="A sequence code from RDKit (currently only supports peptides).",
    )
    fasta: Optional[Union[ChemCode, str]] = Field(
        None, description="A FASTA code (currently only supports peptides)."
    )
    helm: Optional[Union[ChemCode, str]] = Field(
        None, description="A HELM code (currently only supports peptides)."
    )


class Molecule(ProtoModel):
    """A representation of a Molecule in MM. This model contains data for symbols, geometry, connectivity, charges,
    residues, etc. while also supporting a wide array of I/O and manipulation capabilities. Charges, masses, geometry,
    and velocities are truncated to 4, 6, 8, and 8 decimal places, respectively, to assist with duplicate detection.
    """

    schema_name: constr(
        strip_whitespace=True, regex="^(mmschema_molecule)$"
    ) = Field(  # type: ignore
        mmschema_molecule_default,
        description=(
            f"The MMSchema specification to which this model conforms. Explicitly fixed as {mmschema_molecule_default}."
        ),
    )
    schema_version: int = Field(  # type: ignore
        0,
        description="The version number of ``schema_name`` to which this model conforms.",
    )
    symbols: Optional[List[str]] = Field(  # type: ignore
        None,
        description="An ordered (natom,) array-like object of atomic elemental symbols. The index of "
        "this attribute sets atomic order for all other per-atom setting like ``real`` and the first "
        "dimension of ``geometry``. Ghost/Virtual atoms must have an entry in this array-like and are "
        "indicated by the matching the 0-indexed indices in ``real`` field.",
    )
    # Basics
    name: Optional[str] = Field(  # type: ignore
        None,
        description="Common or human-readable name to assign to this molecule. This field can be arbitrary; see "
        "``identifiers`` for well-defined labels.",
    )
    identifiers: Optional[Identifiers] = Field(  # type: ignore
        None,
        description="An optional dictionary of additional identifiers by which this molecule can be referenced, "
        "such as INCHI, canonical SMILES, etc. See the :class:``Identifiers`` model for more details.",
    )
    comment: Optional[str] = Field(  # type: ignore
        None,
        description="Additional comments for this molecule. Intended for pure human/user consumption and clarity.",
    )
    ndim: Optional[int] = Field(  # type: ignore
        3, description="Number of spatial dimensions."
    )
    # Molecular data
    real_: Optional[Array[bool]] = Field(  # type: ignore
        None,
        description="The ordered array indicating if each atom is real (``True``) or "
        "ghost/virtual (``False``). Index matches the 0-indexed indices of all other per-atom settings like "
        "``symbols`` and the first dimension of ``geometry``. If this is not provided, all atoms are assumed "
        "to be real (``True``). If this is provided, the reality or ghostedness of every atom must be specified.",
        shape=["nat"],
    )
    atom_labels_: Optional[List[str]] = Field(  # type: ignore
        None,
        description="Additional per-atom labels as an array of strings. Typical use is in "
        "model conversions, such as Elemental <-> Molpro and not typically something which should be user "
        "assigned. See the ``comments`` field for general human-consumable text to affix to the molecule.",
        shape=["nat"],
    )
    atomic_numbers_: Optional[Array[numpy.int16]] = Field(  # type: ignore
        None,
        description="An optional ordered 1-D array-like object of atomic numbers of shape (nat,). Index "
        "matches the 0-indexed indices of all other per-atom settings like ``symbols`` and ``real``. "
        "Values are inferred from the ``symbols`` list if not explicitly set. "
        "Ghostedness should be indicated through ``real`` field, not zeros here.",
        shape=["nat"],
    )
    mass_numbers_: Optional[Array[numpy.int16]] = Field(  # type: ignore
        None,
        description="An optional ordered 1-D array-like object of atomic *mass* numbers of shape (nat). Index "
        "matches the 0-indexed indices of all other per-atom settings like ``symbols`` and ``real``. "
        "Values are inferred from the most common isotopes of the ``symbols`` list if not explicitly set. "
        "If single isotope not (yet) known for an atom, -1 is placeholder.",
        shape=["nat"],
    )
    masses_: Optional[Array[float]] = Field(  # type: ignore
        None,
        description="The ordered array of particle masses. Index order "
        "matches the 0-indexed indices of all other per-atom fields like ``symbols`` and ``real``. If "
        "this is not provided, the mass of each atom is inferred from its most common isotope. If this "
        "is provided, it must be the same length as ``symbols`` but can accept ``None`` entries for "
        "standard masses to infer from the same index in the ``symbols`` field.",
        shape=["nat"],
    )
    masses_units: Optional[str] = Field(  # type: ignore
        "amu",
        description="Units for atomic masses. Defaults to unified atomic mass unit.",
    )
    molecular_charge: float = Field(  # type: ignore
        0.0,
        description="The net electrostatic charge of the molecule. Default unit is electron Volt.",
    )
    molecular_charge_units: Optional[str] = Field(  # type: ignore
        "e", description="Units for molecular charge. Defaults to electron Volt."
    )
    geometry: Optional[Array[float]] = Field(  # type: ignore
        None,
        description="An ordered (natom*ndim,) array for XYZ atomic coordinates. Default unit is Angstrom.",
    )
    geometry_units: Optional[str] = Field(  # type: ignore
        "angstrom", description="Units for atomic geometry. Defaults to Angstroms."
    )
    velocities: Optional[Array[float]] = Field(  # type: ignore
        None,
        description="An ordered (natoms*ndim,) array for XYZ atomic velocities. Default unit is "
        "Angstroms/femtoseconds.",
    )
    velocities_units: Optional[str] = Field(  # type: ignore
        "angstrom/fs",
        description="Units for atomic velocities. Defaults to Angstroms/femtoseconds.",
    )
    # Topological data
    connectivity_: Optional[List[Tuple[int, int, float]]] = Field(  # type: ignore
        None,
        description="A list of bonds within the molecule. Each entry is a tuple "
        "of ``(atom_index_A, atom_index_B, bond_order)`` where the ``atom_index`` "
        "matches the 0-indexed indices of all other per-atom settings like ``symbols`` and ``real``. "
        "Bonds may be freely reordered and inverted.",
        min_items=1,
    )
    substructs: Optional[List[Tuple[str, int]]] = Field(  # type: ignore
        None,
        description="A list of (name, num) of connected atoms constituting the building block (e.g. monomer) "
        "of the structure (e.g. a polymer). Order follows atomic indices from 0 till Natoms-1. E.g. [('ALA', 4), ...] "
        "means atom1 belongs to aminoacid alanine with residue number 4.",
    )
    chains: Optional[Dict[str, List[int]]] = Field(  # type: ignore
        None,
        description="A sequence of connected substructures forming a subunit that is not bonded to any "
        "other subunit. For example, a hemoglobin molecule consists of four chains that are not connected to one another.",
    )
    # Extras
    provenance: Provenance = Field(
        default_factory=functools.partial(provenance_stamp, __name__),
        description="The provenance information about how this Molecule (and its attributes) were generated, "
        "provided, and manipulated.",
    )
    extras: Dict[str, Any] = Field(  # type: ignore
        None,
        description="Additional information to bundle with the molecule. Use for schema development and scratch space.",
    )

    class Config(ProtoModel.Config):
        serialize_skip_defaults = True
        repr_style = lambda self: [("name", self.name), ("hash", self.get_hash()[:7])]
        fields = {
            "masses_": "masses",
            "real_": "real",
            "atom_labels_": "atom_labels",
            "atomic_numbers_": "atomic_numbers",
            "mass_numbers_": "mass_numbers",
            "connectivity_": "connectivity",
            # below addresses the draft-04 issue until https://github.com/samuelcolvin/pydantic/issues/1478 .
        }
        schema_extra = "http://json-schema.org/draft-04/schema#"

    def __init__(self, **kwargs: Optional[Dict[str, Any]]) -> None:
        """
        Initializes the molecule object from dictionary-like values.
        Parameters
        ----------
        **kwargs : Any
            The values of the Molecule object attributes.
        """
        kwargs["schema_name"] = kwargs.pop("schema_name", "mmschema_molecule")
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

        if values.get("symbols") is None:
            raise ValueError(
                "Either symbols or atomic_numbers must be supplied for a unique definition of a Molecule."
            )

        if not values.get("name"):
            from qcelemental.molparse.to_string import formula_generator

            values["name"] = formula_generator(values["symbols"])

    # Validators
    @validator("*", pre=True)
    def _empty_must_none(cls, v, values):
        """
        Makes sure empty lists or tuples are converted to None.
        """
        if isinstance(v, List) or isinstance(v, Tuple):
            return v if v else None
        return v

    @validator("geometry", "velocities")
    def _must_be_3n(cls, v, values, **kwargs):
        if v is not None:
            n = len(values["symbols"])
            try:
                v.reshape(n, values["ndim"])
            except (ValueError, AttributeError):
                raise ValueError("Array must be castable to shape (natom,ndim)!")
        return v

    @validator("masses_", "real_")
    def _must_be_n(cls, v, values, **kwargs):
        n = len(values["symbols"])
        if len(v) != n:
            raise ValueError(
                "Masses and Real must be same number of entries as Symbols"
            )
        return v

    @validator("real_")
    def _populate_real(cls, v, values, **kwargs):
        # Can't use geometry here since its already been validated and not in values
        n = len(values["symbols"])
        if len(v) == 0:
            v = numpy.array([True for _ in range(n)])
        return v

    @validator("geometry")
    def _valid_dims(cls, v, values, **kwargs):
        n = len(values["symbols"])
        try:
            v = v.reshape(n, values["ndim"])
        except (ValueError, AttributeError):
            raise ValueError(
                f"Geometry must be castable to shape (Natoms,{values['ndim']})!"
            )
        return v

    # Properties
    @property
    def masses(self) -> Array[float]:
        masses = self.__dict__.get("masses_")
        if masses is None:
            try:
                masses = numpy.array(
                    [qcelemental.periodictable.to_mass(x) for x in self.symbols]
                )
            except Exception:
                masses = None
        return masses

    @property
    def real(self) -> Array[bool]:
        real = self.__dict__.get("real_")
        if real is None:
            real = numpy.array([True for x in self.symbols])
        return real

    @property
    def atom_labels(self) -> Array[str]:
        atom_labels = self.__dict__.get("atom_labels_")
        if atom_labels is None:
            atom_labels = numpy.array(["" for x in self.symbols])
        return atom_labels

    @property
    def atomic_numbers(self) -> Array[numpy.int16]:
        atomic_numbers = self.__dict__.get("atomic_numbers_")
        if atomic_numbers is None:
            try:
                atomic_numbers = numpy.array(
                    [qcelemental.periodictable.to_Z(x) for x in self.symbols]
                )
            except Exception:
                atomic_numbers = None
        return atomic_numbers

    @property
    def mass_numbers(self) -> Array[numpy.int16]:
        mass_numbers = self.__dict__.get("mass_numbers_")
        if mass_numbers is None:
            mass_numbers = numpy.array(
                [qcelemental.periodictable.to_A(x) for x in self.symbols]
            )
        return mass_numbers

    @property
    def connectivity(self) -> List[Tuple[int, int, float]]:
        connectivity = self.__dict__.get("connectivity_")
        # default is None, not []
        return connectivity

    @property
    def units(self):
        return {
            val.alias: val.field_info.extra.get("units")
            for key, val in self.__fields__.items()
            if "units" in val.field_info.extra
        }

    @property
    def hash_fields(self):
        return [
            "symbols",
            "masses",
            "molecular_charge",
            "real",
            "geometry",
            "velocities",
            "connectivity",
        ]

    # Representation -> used by qcelemental's __repr__
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

    def pretty_print(self):
        """Prints the molecule. Not sure yet what this is used for, but I have modified the
        original implementation from qcelemental.
        Returns
        -------
        str
            Molecule representation in terms of its total charge and coordinates (geometry).
        """
        text = ""

        text += """    Geometry (in {0:s}), charge = {1:.1f} (in {2:s}):\n\n""".format(
            self.geometry_units, self.molecular_charge, self.molecular_charge_units
        )
        text += """       Center              X                  Y                   Z       \n"""
        text += """    ------------   -----------------  -----------------  -----------------\n"""

        for i in range(len(self.geometry)):
            text += """    {0:8s}{1:4s} """.format(
                self.symbols[i], "" if self.real[i] else "(Gh)"
            )
            for j in range(self.geometry.shape(axis=1)):
                text += """  {0:17.12f}""".format(self.geometry[i][j])
            text += "\n"

        return text

    def get_molecular_formula(self, order: Optional[str] = "alphabetical") -> str:
        """
        Returns the molecular formula for a molecule.

        Parameters
        ----------
        order: str, optional
            Sorting order of the formula. Valid choices are "alphabetical" and "hill".

        Returns
        -------
        str
            The molecular formula.

        Examples
        --------

        >>> methane = qcelemental.models.Molecule('''
        ... H      0.5288      0.1610      0.9359
        ... C      0.0000      0.0000      0.0000
        ... H      0.2051      0.8240     -0.6786
        ... H      0.3345     -0.9314     -0.4496
        ... H     -1.0685     -0.0537      0.1921
        ... ''')
        >>> methane.get_molecular_formula()
        CH4

        >>> hcl = qcelemental.models.Molecule('''
        ... H      0.0000      0.0000      0.0000
        ... Cl     0.0000      0.0000      1.2000
        ... ''')
        >>> hcl.get_molecular_formula()
        ClH

        """
        return qcelemental.molutil.molecular_formula_from_symbols(
            symbols=self.symbols, order=order
        )

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
                if field == "geometry":
                    data = qcelemental.models.molecule.float_prep(data, GEOMETRY_NOISE)
                elif field == "velocities":
                    data = qcelemental.models.molecule.float_prep(data, VELOCITY_NOISE)
                elif field == "fragment_charges":
                    data = qcelemental.models.molecule.float_prep(data, CHARGE_NOISE)
                elif field == "molecular_charge":
                    data = qcelemental.models.molecule.float_prep(data, CHARGE_NOISE)
                elif field == "masses":
                    data = qcelemental.models.molecule.float_prep(data, MASS_NOISE)

                concat += json.dumps(data, default=lambda x: x.ravel().tolist())

        m.update(concat.encode("utf-8"))
        return m.hexdigest()

    # Constructors
    @classmethod
    def from_file(
        cls,
        filename: str,
        top_filename: Optional[str] = None,
        dtype: Optional[str] = None,
        *,
        translator: Optional[str] = None,
        **kwargs: Optional[Dict[str, Any]],
    ) -> "Molecule":
        """
        Constructs a Molecule object from a file.
        Parameters
        ----------
        filename : str
            The molecular structure filename to read
        top_filename: str, optional
            The topology i.e. connectivity filename to read
        dtype : str, optional
            The type of file to interpret. If not set, mmelemental attempts to discover the file type.
        translator: Optional[str], optional
            Translator name e.g. mmic_rdkit. Takes precedence over dtype. If unset,
            MMElemental attempts to find an appropriate translator if it is registered
            in the :class:``TransComponent`` class.
        **kwargs: Optional[Dict[str, Any]], optional
            Any additional keywords to pass to the constructor
        Returns
        -------
        Molecule
            A constructed Molecule class.
        """
        file_ext = Path(filename).suffix if filename else None

        if file_ext in qcelemental.models.molecule._extension_map:
            if top_filename:
                raise TypeError(
                    "Molecule topology must be supplied in a single JSON (or similar) file."
                )

            if dtype is None:
                dtype = qcelemental.models.molecule._extension_map[file_ext]

            # Raw string type, read and pass through
            if dtype == "json":
                with open(filename, "r") as infile:
                    data = json.load(infile)
                dtype = "dict"
            else:
                raise KeyError(f"Data type not supported: {dtype}.")

            return cls.from_data(data, dtype=dtype, **kwargs)

        fileobj = FileOutput(path=filename) if filename else None
        top_fileobj = FileOutput(path=top_filename) if top_filename else None

        dtype = dtype or fileobj.ext.strip(".")

        # Generic translator component
        try:
            from mmic_translator.components import TransComponent
        except Exception:
            TransComponent = None

        if not translator:
            if not TransComponent:
                raise ModuleNotFoundError(_trans_nfound_msg)
            inst_trans = TransComponent.installed_comps()

            while not translator:
                translator = TransComponent.find_molread_tk(
                    fileobj.ext, trans=inst_trans
                )
                if not translator:
                    raise ValueError(
                        f"There is no installed translator for reading file {filename}. Please install an appropriate translator."
                    )

                # We should be able to always import the translator
                mod = importlib.import_module(translator)

                # If top if supplied, make sure the translator supports the top file extension
                if top_fileobj:
                    top_ext = top_fileobj.ext
                    if top_ext not in mod.ffread_ext_maps:
                        inst_trans.remove(translator)
                        translator = None
                    if len(inst_trans) == 0:
                        raise ValueError(
                            f"There is no installed translator for concurrently reading files {filename} and {top_filename}.\n"
                            + "Please install an appropriate translator."
                        )
        elif importlib.util.find_spec(translator):
            mod = importlib.import_module(translator)

        tkmol_class = mod._classes_map.get("Molecule")

        if not tkmol_class:
            raise ValueError(
                f"No Molecule model found while looking in translator: {translator}."
            )

        tkmol = tkmol_class.from_file(
            filename=fileobj.abs_path if fileobj else None,
            top_filename=top_fileobj.abs_path if top_fileobj else None,
            dtype=dtype,
        )

        return cls.from_data(tkmol, dtype=tkmol.dtype, **kwargs)

    @classmethod
    def from_data(
        cls,
        data: Optional[Any] = None,
        dtype: Optional[str] = None,
        **kwargs: Dict[str, Any],
    ) -> "Molecule":
        """
        Constructs a Molecule object from a data object.
        Parameters
        ----------
        data: Any, optional
            Data to construct Molecule from such as a data object (e.g. MDAnalysis.Universe) or dict.
        dtype: str, optional
            How to interpret the data, if not passed attempts to discover this based on input type.
        **kwargs: Optional[Dict[str, Any]], optional
            Additional kwargs to pass to the constructors.
        Returns
        -------
        Molecule
            A constructed Molecule class.
        """
        if isinstance(data, str):
            if not dtype:
                raise ValueError(
                    "You must supply dtype for proper interpretation of symbolic data e.g. MDAnalysis, smiles, etc."
                )
            code = ChemCode(code=data, dtype=dtype)
            symbols = list(code.code)  # this is garbage, must be replaced
            return Molecule(identifiers={dtype: data}, symbols=symbols)

            try:
                import mmic_molconv
            except Exception as e:
                raise ValueError(f"Failed in importing mmic_molconv. Exception: {e}")

            return mmic_molconv.RunComponent.compute(
                {"data": data, "kwargs": kwargs}
            ).molecule

        elif isinstance(data, dict):
            kwargs.pop("dtype", None)  # remove dtype if supplied
            kwargs.update(data)
            return cls(**kwargs)

        return data.to_schema(**kwargs)

    def to_file(
        self,
        filename: str,
        dtype: Optional[str] = None,
        *,
        translator: Optional[str] = None,
        **kwargs: Optional[Dict[str, Any]],
    ) -> None:
        """Writes the Molecule to a file.
        Parameters
        ----------
        filename : str
            The filename to write to
        dtype : Optional[str], optional
            The type of file to write (e.g. json, pdb, etc.), attempts to infer dtype from
            file extension if not provided.
        translator: Optional[str], optional
            Translator name e.g. mmic_rdkit. Takes precedence over dtype. If unset,
            MMElemental attempts to find an appropriate translator if it is registered
            in the :class:``TransComponent`` class.
        **kwargs: Optional[Dict[str, Any]], optional
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
        else:  # look for an installed mmic_translator
            try:
                from mmic_translator.components import TransComponent
            except Exception:
                TransComponent = None

            if not TransComponent:
                raise ModuleNotFoundError(_trans_nfound_msg)
            translator = TransComponent.find_molwrite_tk(ext)

            if not translator:
                raise NotImplementedError(
                    f"File extension {ext} not supported with any installed translators."
                )

            tkmol = self.to_data(translator=translator, **kwargs)
            tkmol.to_file(filename, dtype=dtype, **kwargs)  # pass dtype?

    def to_data(
        self,
        dtype: Optional[str] = None,
        *,
        translator: Optional[str] = None,
        **kwargs: Optional[Dict[str, Any]],
    ) -> "ToolkitModel":
        """Converts Molecule to toolkit-specific molecule (e.g. rdkit, MDAnalysis, parmed).
        Parameters
        ----------
        dtype: Optional[str], optional
            The type of data object to convert to e.g. MDAnalysis, rdkit, parmed, etc.
        translator: Optional[str], optional
            Translator name e.g. mmic_rdkit. Takes precedence over dtype. If unset,
            MMElemental attempts to find an appropriate translator if it is registered
            in the :class:``TransComponent`` class.
        **kwargs: Optional[Dict[str, Any]], optional
            Additional kwargs to pass to the constructor.
        Returns
        -------
        ToolkitModel
            Toolkit-specific molecule model
        """
        try:
            from mmic_translator.components import TransComponent
        except Exception:
            TransComponent = None

        if not translator:
            if not TransComponent:
                raise ModuleNotFoundError(_trans_nfound_msg)
            if not dtype:
                raise ValueError(
                    f"Either translator or dtype must be supplied when calling {__name__}."
                )
            translator = TransComponent.find_trans(dtype)

        if importlib.util.find_spec(translator):
            mod = importlib.import_module(translator)
            tkmol = mod._classes_map.get("Molecule")

            if not tkmol:
                raise ValueError(
                    f"No Molecule model found while looking in translator: {translator}."
                )

            return tkmol.from_schema(self)
        else:
            raise NotImplementedError(
                f"translator {translator} not available. Make sure it is properly installed."
            )
