import numpy
from typing import List, Optional, Any, Dict, Union
from pydantic import Field, constr, validator, root_validator
import importlib
from pathlib import Path
import hashlib
import json
from functools import partial

# MM models
from mmelemental.models.util.output import FileOutput
from .topology import Topology
from mmelemental.models.base import ProtoModel, Provenance, provenance_stamp
from cmselemental.types import Array
from mmelemental.util.data import (
    float_prep,
    NUMPY_UNI,
    NUMPY_INT,
    NUMPY_FLOAT,
    GEOMETRY_NOISE,
    VELOCITY_NOISE,
    MASS_NOISE,
    CHARGE_NOISE,
)
from mmelemental.util.units import (
    LENGTH_DIM,
    MASS_DIM,
    TIME_DIM,
    CURRENT_DIM,
)
from cmselemental.util import yaml_import, which_import

__all__ = ["Molecule"]


_trans_nfound_msg = "MMElemental translation requires mmic & mmic_translator. \
Solve by: pip install mmic mmic_translator"
_qcel_nfound_msg = "MMElemental feature requires qcelemental. \
Solve by: pip install qcelemental"
mmschema_molecule_default = "mmschema_molecule"


class Identifiers(ProtoModel):
    """
    An extension of the qcelemental.models.molecule.Identifiers for RDKit constructors.
    See `link <https://rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html>`_ for more info.
    """

    molecule_hash: Optional[str] = None
    molecular_formula: Optional[str] = None
    smiles: Optional[str] = None
    inchi: Optional[str] = None
    inchikey: Optional[str] = None
    canonical_explicit_hydrogen_smiles: Optional[str] = None
    canonical_isomeric_explicit_hydrogen_mapped_smiles: Optional[str] = None
    canonical_isomeric_explicit_hydrogen_smiles: Optional[str] = None
    canonical_isomeric_smiles: Optional[str] = None
    canonical_smiles: Optional[str] = None
    pubchem_cid: Optional[str] = Field(None, description="PubChem Compound ID")
    pubchem_sid: Optional[str] = Field(None, description="PubChem Substance ID")
    pubchem_conformerid: Optional[str] = Field(None, description="PubChem Conformer ID")
    smiles: Optional[str] = Field(
        None, description="A simplified molecular-input line-entry system code."
    )
    smarts: Optional[str] = Field(
        None,
        description="A SMILES arbitrary target specification code for defining substructures.",
    )
    inchi: Optional[str] = Field(
        None, description="An international chemical identifier code."
    )
    sequence: Optional[str] = Field(
        None,
        description="A sequence code from RDKit (currently only supports peptides).",
    )
    fasta: Optional[str] = Field(
        None, description="A FASTA code (currently only supports peptides)."
    )
    helm: Optional[str] = Field(
        None, description="A HELM code (currently only supports peptides)."
    )

    class Config(ProtoModel.Config):
        serialize_skip_defaults = True


class Molecule(ProtoModel):
    """A representation of a Molecule in MM. This model contains data for symbols, geometry, connectivity, charges,
    residues, etc. while also supporting a wide array of I/O and manipulation capabilities. Charges, masses, geometry,
    and velocities are truncated to 4, 6, 8, and 8 decimal places, respectively, to assist with duplicate detection.
    """

    # Basic fields
    schema_name: constr(
        strip_whitespace=True, regex=mmschema_molecule_default
    ) = Field(  # type: ignore
        mmschema_molecule_default,
        description=(
            f"The MMSchema specification to which this model conforms. Explicitly fixed as {mmschema_molecule_default}."
        ),
    )
    schema_version: int = Field(  # type: ignore
        1,
        description="The version number of ``schema_name`` to which this model conforms.",
    )
    symbols: Optional[Array[str]] = Field(  # type: ignore
        None,
        description="An ordered (natom,) array-like object of particle symbols. The index of "
        "this attribute sets the order for all other per-particle setting like ``geometry`` and the first "
        "dimension of ``geometry``.",
    )
    name: Optional[str] = Field(  # type: ignore
        None,
        description="Common or human-readable name to assign to this molecule. This field can be arbitrary.",
    )
    identifiers: Optional[Identifiers] = Field(  # type: ignore
        None,
        description="An optional dictionary of additional identifiers by which this molecule can be referenced, "
        "such as INCHI, canonical SMILES, etc. See the :class:`Identifiers` model for more details.",
    )
    comment: Optional[str] = Field(  # type: ignore
        None,
        description="Additional comments for this molecule. Intended for pure human/user consumption and clarity.",
    )
    ndim: Optional[int] = Field(  # type: ignore
        3, description="Number of spatial dimensions."
    )
    # Molecular data
    atom_labels: Optional[Array[str]] = Field(  # type: ignore
        None,
        description="Additional per-atom labels as an array of strings. Typical use is in "
        "model conversions, such as Elemental <-> Molpro and not typically something which should be user "
        "assigned. See the ``comments`` field for general human-consumable text to affix to the molecule.",
    )
    atomic_numbers_: Optional[Array[numpy.dtype(NUMPY_INT)]] = Field(  # type: ignore
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
    masses_: Optional[Array[numpy.dtype(NUMPY_FLOAT)]] = Field(  # type: ignore
        None,
        description="The ordered array of particle masses. Index order "
        "matches the 0-indexed indices of all other per-atom fields like ``symbols`` and ``real``. If "
        "this is not provided, the mass of each atom is inferred from its most common isotope. If this "
        "is provided, it must be the same length as ``symbols`` but can accept ``None`` entries for "
        "standard masses to infer from the same index in the ``symbols`` field.",
    )
    masses_units: Optional[str] = Field(  # type: ignore
        "unified_atomic_mass_unit",
        description="Units for atomic masses. Defaults to unified atomic mass unit.",
        dimensionality=MASS_DIM,
    )
    molecular_charge: Optional[float] = Field(  # type: ignore
        0.0,
        description="The net electrostatic charge of the molecule. Default unit is elementary charge.",
    )
    molecular_charge_units: Optional[str] = Field(  # type: ignore
        "elementary_charge",
        description="Units for molecular charge. Defaults to elementary charge.",
        dimensionality=CURRENT_DIM * TIME_DIM,
    )
    formal_charges: Optional[Array[numpy.dtype(NUMPY_FLOAT)]] = Field(  # type: ignore
        None, description="Formal charges of all particles/atoms."
    )
    formal_charges_units: Optional[str] = Field(  # type: ignore
        "elementary_charge",
        description="Units for formal charges. Defaults to elementary charge.",
        dimensionality=CURRENT_DIM * TIME_DIM,
    )
    partial_charges: Optional[Array[numpy.dtype(NUMPY_FLOAT)]] = Field(  # type: ignore
        None, description="Assigned partial charges of all particles/atoms."
    )
    partial_charges_units: Optional[str] = Field(  # type: ignore
        "elementary_charge",
        description="Units for partial charges. Defaults to elementary charge.",
        dimensionality=CURRENT_DIM * TIME_DIM,
    )
    geometry: Optional[Array[numpy.dtype(NUMPY_FLOAT)]] = Field(  # type: ignore
        None,
        description="An ordered (natom*ndim,) array for XYZ atomic coordinates. Default unit is Angstrom.",
    )
    geometry_units: Optional[str] = Field(  # type: ignore
        "angstrom",
        description="Units for atomic geometry. Defaults to Angstroms.",
        dimensionality=LENGTH_DIM,
    )
    velocities: Optional[Array[NUMPY_FLOAT]] = Field(  # type: ignore
        None,
        description="An ordered (natoms*ndim,) array for XYZ atomic velocities. Default unit is "
        "angstroms/femtoseconds.",
    )
    velocities_units: Optional[str] = Field(  # type: ignore
        "angstrom / femtosecond",
        description="Units for atomic velocities. Defaults to Angstroms/femtoseconds.",
        dimensionality=LENGTH_DIM / TIME_DIM,
    )
    # Topological data
    connectivity: Optional[Array[numpy.dtype(f"{NUMPY_INT}, {NUMPY_INT}, {NUMPY_FLOAT}")]] = Field(  # type: ignore
        None,
        description="A list of bonds within the molecule. Each entry is a tuple "
        "of ``(atom_index_A, atom_index_B, bond_order)`` where the ``atom_index`` "
        "matches the 0-indexed indices of all other per-atom settings like ``symbols`` and ``real``. "
        "Bonds may be freely reordered and inverted.",
    )
    substructs: Optional[Array[numpy.dtype(f"{NUMPY_UNI}, {NUMPY_INT}")]] = Field(  # type: ignore
        None,
        description="A list of (name, num) of connected atoms constituting the building block (e.g. monomer) "
        "of the structure (e.g. a polymer). Order follows atomic indices from 0 till Natoms-1. E.g. [('ALA', 4), ...] "
        "means atom1 belongs to aminoacid alanine with residue number 4. Substruct name is max 4 characters.",
    )
    # Extras
    provenance: Provenance = Field(
        default_factory=partial(provenance_stamp, __name__),
        description="The provenance information about how this object (and its attributes) were generated, "
        "provided, and manipulated.",
    )
    extras: Optional[Dict[str, Any]] = Field(  # type: ignore
        None,
        description="Additional information to bundle with the molecule. Use for schema development and scratch space.",
    )
    hash: str = Field(
        None,
        description="The hash code that unique identifies this object. Typically not manually assigned but left for "
        "MMElemental to handle.",
    )

    class Config(ProtoModel.Config):
        repr_style = lambda self: [("name", self.name), ("hash", self.get_hash()[:7])]
        fields = {
            "masses_": "masses",
            "atomic_numbers_": "atomic_numbers",
        }

    def __init__(self, **kwargs: Optional[Dict[str, Any]]) -> None:
        """
        Initializes the molecule object from dictionary-like values.

        Parameters
        ----------
        **kwargs : Any
            The values of the Molecule object attributes.

        """
        atomic_numbers = kwargs.get("atomic_numbers")
        if atomic_numbers is not None:
            if not which_import("qcelemental", return_bool=True):
                raise ModuleNotFoundError(_qcel_nfound_msg)

            import qcelemental

            if kwargs.get("symbols") is None:

                kwargs["symbols"] = [
                    qcelemental.periodictable.to_E(x) for x in atomic_numbers
                ]

        # We are pulling out the values *explicitly* so that the pydantic skip_defaults works as expected
        # All attributes set below are equivalent to the default set.
        super().__init__(**kwargs)

        values = self.__dict__

        if not values.get("name"):
            if not which_import("qcelemental", return_bool=True):
                raise ModuleNotFoundError(_qcel_nfound_msg)

            from qcelemental.molparse.to_string import formula_generator

            values["name"] = formula_generator(values["symbols"])

        if values.get("hash") is None:
            values["hash"] = self.get_hash()
        else:
            assert (
                values["hash"] == self.get_hash()
            ), "Model data inconsistent with stored hash code!"

        if values.get("symbols") is None:
            raise ValueError(
                "Either symbols or atomic_numbers must be supplied for a unique definition of a Molecule."
            )

    # Validators
    @validator("*", pre=True)
    def _empty_must_none(cls, v, values):
        """
        Makes sure empty lists, tuples, or arrays are converted to None.
        """
        if isinstance(v, (list, tuple, numpy.ndarray)):
            return v if len(v) else None
        return v

    @validator("substructs", "connectivity", pre=True)
    def _valid_tuple_array(cls, v):
        assert isinstance(
            v, (list, tuple, numpy.ndarray)
        ), "Substructs must be a list, tuple, or array."
        if isinstance(v[0], list):
            return [tuple(item) for item in v]
        return v

    @validator("geometry", "velocities")
    def _valid_shape(cls, v, values, **kwargs):
        if v is not None:
            n = len(values["symbols"])
            try:
                v.reshape(n, values["ndim"])
            except (ValueError, AttributeError):
                raise ValueError("Array must be castable to shape (natom,ndim)!")
            if v.ndim > 1:
                assert (
                    v.ndim == 2
                ), f"Array must be either 1D or 2D! Given ndim: {v.ndim}."
                assert (
                    v.shape[1] == values["ndim"]
                ), f"Array column number ({v.shape[1]}) should be equal to ndim ({values['ndim']})."
            return v.flatten()

    @root_validator
    def _must_be_n(cls, values):
        assert values["symbols"] is not None, "symbols is required."
        n = len(values["symbols"])
        for key in [
            "masses",
            "substructs",
            "atom_labels",
            "atomic_numbers",
            "mass_numbers",
        ]:
            if values.get(key) is not None:
                assert (
                    len(values[key]) == n
                ), f"{key} (len(values[key])) must have same number of entries (n) as Symbols."
        return values

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
    def atomic_numbers(self) -> Array[numpy.int16]:
        atomic_numbers = self.__dict__.get("atomic_numbers_")
        if atomic_numbers is None:
            try:
                import qcelemental

                atomic_numbers = numpy.array(
                    [qcelemental.periodictable.to_Z(x) for x in self.symbols]
                )
            except Exception:
                atomic_numbers = None
        return atomic_numbers

    @property
    def hash_fields(self):
        return [
            "symbols",
            "masses",
            "masses_units",
            "molecular_charge",
            "molecular_charge_units",
            "geometry",
            "geometry_units",
            "velocities",
            "velocities_units",
            "connectivity",
        ]

    def _ipython_display_(self, **kwargs) -> None:
        try:
            self.show()._ipython_display_(**kwargs)
        except ModuleNotFoundError:
            from IPython.display import display

            display(f"Install nglview for interactive visualization.", f"{repr(self)}")

    def __repr_args__(self) -> "ReprArgs":
        return [("name", self.name), ("hash", self.get_hash()[:7])]

    def __hash__(self):
        return hash(self.get_hash())

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

    def get_topology(self):

        return Topology(
            name="top_from_mol",
            schema_version=self.schema_version,
            symbols=self.symbols,
            atom_labels=self.atom_labels,
            atomic_numbers=self.atomic_numbers,
            mass_numbers=self.mass_numbers,
            masses=self.masses,
            masses_units=self.masses_units,
            molecular_charge=self.molecular_charge,
            molecular_charge_units=self.molecular_charge_units,
            connectivity=self.connectivity,
        )

    def get_substructs(self) -> List:
        """Removes duplicate entries in substructs while preserving the order."""
        seen = set()
        seen_add = seen.add
        return [x for x in self.substructs.tolist() if not (x in seen or seen_add(x))]

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

        """
        if not which_import("qcelemental", return_bool=True):
            raise ModuleNotFoundError(_qcel_nfound_msg)
        import qcelemental

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
                    data = float_prep(data, GEOMETRY_NOISE)
                elif field == "velocities":
                    data = float_prep(data, VELOCITY_NOISE)
                elif field == "molecular_charge":
                    data = float_prep(data, CHARGE_NOISE)
                elif field == "masses":
                    data = float_prep(data, MASS_NOISE)

                concat += json.dumps(
                    data, default=lambda x: x.ravel().tolist()
                )  # if serialization fails, assume type is numpy.ndarray

        m.update(concat.encode("utf-8"))
        return m.hexdigest()

    def show(self, ngl_kwargs: Optional[Dict[str, Any]] = None) -> "nglview.NGLWidget":  # type: ignore
        r"""Creates a 3D representation of a moleucle that can be manipulated in Jupyter Notebooks and exported as
        images (`.png`).

        Parameters
        ----------
        ngl_kwargs
            Addition nglview NGLWidget kwargs.

        Returns
        -------
        nglview.NGLWidget
            An nglview view of the molecule.

        """

        if not which_import("nglview", return_bool=True):
            raise ModuleNotFoundError(
                f"Python module nglwview not found. Solve by installing it: `pip install nglview`"
            )  # pragma: no cover

        import nglview  # type: ignore

        if ngl_kwargs is None:
            ngl_kwargs = {}

        structure = nglview.adaptor.MMElementalStructure(self)
        widget = nglview.NGLWidget(structure, **ngl_kwargs)
        return widget

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
            in the :class:`TransComponent` class.
        **kwargs: Optional[Dict[str, Any]], optional
            Any additional keywords to pass to the constructor

        Returns
        -------
        Molecule
            A Molecule object.

        """
        file_ext = Path(filename).suffix if filename else None

        if file_ext in [".json", ".js", ".yaml", ".yml"]:
            if top_filename:
                raise TypeError(
                    "Molecule topology must be supplied in a single JSON (or YAML) file."
                )

            if dtype is None:
                if file_ext in [".json", ".js"]:
                    dtype = "json"
                elif file_ext in [".yaml", ".yml"]:
                    dtype = "yaml"

            # Raw string type, read and pass through
            if dtype == "json":
                with open(filename, "r") as infile:
                    data = json.load(infile)
                dtype = "dict"
            elif dtype == "yaml":
                with open(filename, "r") as infile:
                    yaml = yaml_import(raise_error=True)
                    data = yaml.safe_load(infile)
                dtype = "dict"
            else:
                raise KeyError(f"Data type not supported: {dtype}.")

            return cls.from_data(data, dtype=dtype, **kwargs)

        fileobj = FileOutput(path=filename) if filename else None
        top_fileobj = FileOutput(path=top_filename) if top_filename else None

        dtype = dtype or fileobj.ext.removeprefix(".")

        # Generic translator component
        if not which_import("mmic_translator", return_bool=True):
            raise ModuleNotFoundError(_trans_nfound_msg)

        from mmic_translator.components import TransComponent

        if not translator:
            inst_trans = TransComponent.installed_comps()

            while not translator:
                translator = TransComponent.find_molread_tk(
                    fileobj.ext, trans=inst_trans
                )
                if not translator:
                    if not top_filename:
                        raise ValueError(
                            f"There is no installed translator for reading file {filename}. Please install an appropriate translator."
                        )
                    else:
                        raise ValueError(
                            f"There is no installed translator for concurrently reading files {filename} and {top_filename}.\n"
                            + "Please install an appropriate translator."
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
        elif which_import(translator, return_bool=True):
            mod = importlib.import_module(translator)
        else:
            raise ModuleNotFoundError(f"Translator {translator} is not installed.")

        try:
            tkmol_class = mod._classes_map.get("Molecule")
        except AttributeError as e:
            raise AttributeError(
                f"{translator} is not a valid Molecule translator. Exception raised: "
                + repr(e)
            )

        if kwargs.pop("debug", None):
            kwargs.update(
                extras={
                    "mmic_translator": {
                        "routine": f"{__name__}.{cls.__name__}.{cls.from_file.__name__}",
                        "translator": (translator, mod.__version__),
                        "engine": tkmol_class.engine(),
                        "model": tkmol_class,
                    },
                }
            )

        tkmol = tkmol_class.from_file(
            filename=fileobj.abs_path if fileobj else None,
            top_filename=top_fileobj.abs_path if top_fileobj else None,
            dtype=dtype,
        )

        dtype = TransComponent.get_dtype(translator)
        return cls.from_data(tkmol.data, dtype=dtype, **kwargs)

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
            A Molecule object.

        Examples
        --------
        >>> mol = mmelemental.models.Molecule.from_data(universe, dtype="mdanalysis")
        >>> mol = mmelemental.models.Molecule.from_data(struct, dtype="parmed")

        """
        if isinstance(data, str):

            if not dtype:
                raise ValueError(
                    "You must supply dtype for proper interpretation of string data e.g. smiles, yaml, json, etc."
                )
            elif dtype == "json":
                assert isinstance(data, str)
                input_dict = json.loads(data)
            elif dtype == "yaml":
                yaml = yaml_import(raise_error=True)
                assert isinstance(data, str)
                input_dict = yaml.safe_load(data)
            elif dtype == "smiles":
                input_dict = {
                    "identifiers": {dtype: data},
                    "symbols": list(data),
                }  # move this to mmic_translator? Must be generalized
            else:
                raise NotImplementedError(f"Data type {dtype} not understood.")
            kwargs.pop("dtype", None)
            kwargs.update(input_dict)
            return cls(**kwargs)

        elif isinstance(data, dict):
            kwargs.pop("dtype", None)  # remove dtype if supplied
            kwargs.update(data)
            return cls(**kwargs)

        if not dtype:
            raise ValueError(
                "You must supply dtype for proper interpretation of data objects e.g. mdanalysis, qcschema, etc."
            )

        if not which_import("mmic_translator", return_bool=True):
            raise ModuleNotFoundError(_trans_nfound_msg)

        from mmic_translator.components import TransComponent

        translator = TransComponent.find_trans(dtype)

        if not translator:
            raise NotImplementedError(
                f"Data type {dtype} not supported with any installed translators."
            )

        if which_import(translator, return_bool=True):
            mod = importlib.import_module(translator)
            tkmol_class = mod._classes_map.get("Molecule")
            tkmol = tkmol_class(data=data)

            if not tkmol:
                raise ValueError(
                    f"No Molecule model found while looking in translator: {translator}."
                )

            if kwargs.pop("debug", None):
                kwargs.update(
                    extras={
                        "debug": {
                            "routine": f"{__name__}.{cls.__name__}.{cls.from_file.__name__}",
                            "translator": (translator, mod.__version__),
                            "engine": tkmol_class.engine(),
                            "model": tkmol_class,
                        }
                    }
                )

            return tkmol.to_schema(**kwargs)
        else:
            raise NotImplementedError(
                f"Translator {translator} not available. Make sure it is properly installed."
            )

    def to_file(
        self,
        filename: str,
        dtype: Optional[str] = None,
        *,
        translator: Optional[str] = None,
        **kwargs: Optional[Dict[str, Any]],
    ) -> None:
        """
        Writes the Molecule to a file.

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
            in the :class:`TransComponent` class.
        **kwargs: Optional[Dict[str, Any]], optional
            Additional kwargs to pass to the constructor.

        """
        if not dtype:
            dtype = Path(filename).suffix[1:]

        if dtype in ["json", "js", "yaml", "yml", "hdf5", "h5"]:
            self.write_file(
                filename,
                encoding=dtype,
                **kwargs,
            )
        else:  # look for an installed mmic_translator
            if not translator:
                if not which_import("mmic_translator", return_bool=True):
                    raise ModuleNotFoundError(_trans_nfound_msg)

                from mmic_translator.components import TransComponent

                translator = TransComponent.find_molwrite_tk("." + dtype)

                if not translator:
                    raise NotImplementedError(
                        f"File extension {dtype} not supported with any installed translators."
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
        """
        Converts Molecule to toolkit-specific molecule (e.g. rdkit, MDAnalysis, parmed).

        Parameters
        ----------
        dtype: str, optional
            The type of data object to convert to e.g. mdanalysis, rdkit, parmed, etc.
        translator: Optional[str], optional
            Translator name e.g. mmic_rdkit. Takes precedence over dtype. If unset,
            MMElemental attempts to find an appropriate translator if it is registered
            in the :class:`TransComponent` class.
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

        if which_import(translator, return_bool=True):
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
