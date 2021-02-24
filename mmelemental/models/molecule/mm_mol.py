import qcelemental
from typing import List, Tuple, Optional, Any, Dict, Union
from pydantic import Field, constr, validator
import importlib
from pathlib import Path
import hashlib
import json

# Import MM models
from mmelemental.models.util.output import FileOutput
from mmelemental.models.base import ToolkitModel
from mmelemental.models.chem.codes import ChemCode
from mmelemental.models.base import Provenance, provenance_stamp

# Import MM components
from mmelemental.components.trans import TransComponent


__all__ = ["Molecule"]

# Rounding quantities for hashing
GEOMETRY_NOISE = 8
VELOCITY_NOISE = 8
FORCE_NOISE = 8
MASS_NOISE = 6
CHARGE_NOISE = 4


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


mmschema_molecule_default = "mmschema_molecule"


class Molecule(qcelemental.models.Molecule):
    """A representation of a Molecule in MM based on QCSchema. This model contains data for symbols, geometry,
    connectivity, charges, residues, etc. while also supporting a wide array of I/O and manipulation capabilities.
    Charges, masses, geometry, velocities, and forces are truncated to 4, 6, 8, 8, and 8 decimal places respectively
    to assist with duplicate detection.
    """

    schema_name: constr(
        strip_whitespace=True, regex="^(mmschema_molecule)$"
    ) = Field(  # type: ignore
        mmschema_molecule_default,
        description=(
            f"The MMSchema specification to which this model conforms. Explicitly fixed as {mmschema_molecule_default}."
        ),
    )
    symbols: Optional[qcelemental.models.types.Array[str]] = Field(  # type: ignore
        None,
        description="An ordered (natom,) array-like object of atomic elemental symbols. The index of "
        "this attribute sets atomic order for all other per-atom setting like ``real`` and the first "
        "dimension of ``geometry``. Ghost/Virtual atoms must have an entry in this array-like and are "
        "indicated by the matching the 0-indexed indices in ``real`` field.",
    )
    masses_units: Optional[str] = Field(  # type: ignore
        "amu",
        description="Units for atomic masses. Defaults to unified atomic mass unit.",
    )
    molecular_charge_units: Optional[str] = Field(  # type: ignore
        "eV", description="Units for molecular charge. Defaults to electron Volt."
    )
    geometry: Optional[qcelemental.models.types.Array[float]] = Field(  # type: ignore
        None,
        description="An ordered (natom,3) array-like for XYZ atomic positions in Angstrom. "
        "Can also accept arrays which can be mapped to (natom,3) such as a 1-D list of length 3*natom, "
        "or the serialized version of the array in (3*natom,) shape; all forms will be reshaped to "
        "(natom,3) for this attribute. Default unit is Angstroms.",
    )
    geometry_units: Optional[str] = Field(  # type: ignore
        "angstrom", description="Units for atomic geometry. Defaults to Angstroms."
    )
    velocities: Optional[qcelemental.models.types.Array[float]] = Field(  # type: ignore
        None,
        description="An ordered (natoms,3) array-like for XYZ atomic velocities in Angstrom/ps. "
        "Can also accept arrays which can be mapped to (natoms,3) such as a 1-D list of length 3*natoms, "
        "or the serialized version of the array in (3*natoms,) shape; all forms will be reshaped to "
        "(natoms,3) for this attribute. Default unit is Angstroms/femtoseconds.",
    )
    velocities_units: Optional[str] = Field(  # type: ignore
        "angstrom/fs",
        description="Units for atomic velocities. Defaults to Angstroms/femtoseconds.",
    )
    forces: Optional[qcelemental.models.types.Array[float]] = Field(  # type: ignore
        None,
        description="An ordered (natoms,3) array-like for XYZ atomic velocities in kJ/mol*Angstrom. "
        "Can also accept arrays which can be mapped to (natoms,3) such as a 1-D list of length 3*natoms, "
        "or the serialized version of the array in (3*natoms,) shape; all forms will be reshaped to "
        "(natoms,3) for this attribute. Default unit is KiloJoules/mol.Angstroms.",
    )
    forces_units: Optional[str] = Field(  # type: ignore
        "kJ/(mol*angstrom)",
        description="Units for atomic forces. Defaults to KiloJoules/mol.Angstroms",
    )
    angles: Optional[List[Tuple[int, int, int]]] = Field(
        None, description="Bond angles in degrees for three connected atoms."
    )
    dihedrals: Optional[List[Tuple[int, int, int, int, int]]] = Field(  # type: ignore
        None,
        description="Dihedral/torsion angles in degrees between planes through two sets of three atoms, having two atoms in common.",
    )
    improper_dihedrals: Optional[
        List[Tuple[int, int, int, int, int]]
    ] = Field(  # type: ignore
        None,
        description="Improper dihedral/torsion angles in degrees between planes through two sets of three atoms, having two atoms in common.",
    )
    residues: Optional[List[Tuple[str, int]]] = Field(  # type: ignore
        None,
        description="A list of (residue_name, residue_num) of connected atoms constituting the building block (monomer) "
        "of a polymer. Order follows atomic indices from 0 till Natoms-1. Residue number starts from 1."
        "\n"
        "E.g. ('ALA', 1) means atom 0 belongs to aminoacid alanine with residue number 1.",
    )
    chains: Optional[Dict[str, List[int]]] = Field(  # type: ignore
        None,
        description="A sequence of connected residues (i.e. polymers) forming a subunit that is not bonded to any "
        "other subunit. For example, a hemoglobin molecule consists of four chains that are not connected to one another.",
    )
    segments: Optional[Dict[str, List[int]]] = Field(
        None, description="..."
    )  # type: ignore
    names: Optional[List[str]] = Field(  # type: ignore
        None, description="A list of atomic label names."
    )
    identifiers: Optional[Identifiers] = Field(  # type: ignore
        None,
        description="An optional dictionary of additional identifiers by which this Molecule can be referenced, "
        "such as INCHI, SMILES, SMARTS, etc. See the :class:``Identifiers`` model for more details.",
    )
    connectivity_: Optional[List[Tuple[int, int, float]]] = Field(  # type: ignore
        None,
        description="A list of bonds within the molecule. Each entry is a tuple "
        "of ``(atom_index_A, atom_index_B, bond_order)`` where the ``atom_index`` "
        "matches the 0-indexed indices of all other per-atom settings like ``symbols`` and ``real``. "
        "Bonds may be freely reordered and inverted.",
        min_items=1,
    )
    provenance: Provenance = Field(  # type: ignore
        provenance_stamp(__name__),
        description="The provenance information about how this object (and its attributes) were generated, "
        "provided, and manipulated.",
    )

    def __init__(self, **kwargs: Any) -> None:
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
                from qcelemental import periodic_table

                kwargs["symbols"] = [
                    periodic_table.periodictable.to_E(x) for x in atomic_numbers
                ]

        # We are pulling out the values *explicitly* so that the pydantic skip_defaults works as expected
        # All attributes set below are equivalent to the default set.
        # We must always prevent QCElemental.Molecule from doing any validation, so validate=False.
        super().__init__(orient=False, validate=False, **kwargs)

        values = self.__dict__

        if values.get("symbols") is not None:
            import numpy

            values["symbols"] = numpy.core.defchararray.title(
                self.symbols
            )  # Title case for consistency
        else:
            raise ValueError(
                "Either symbols or atomic_numbers must be supplied for a unique definition of a Molecule."
            )

        if not values.get("name"):
            from qcelemental.molparse.to_string import formula_generator

            values["name"] = formula_generator(values["symbols"])

    @validator("*", pre=True)
    def emptyIsNone(cls, v, values):
        """ 
        Makes sure empty lists or tuples are converted to None.
        """
        if isinstance(v, List) or isinstance(v, Tuple):
            return v if v else None
        return v

    def __eq__(self, other):
        """
        Checks if two molecules are identical. This is a molecular identity defined
        by scientific terms, and not programing terms, so it's less rigorous than
        a programmatic equality or a memory equivalent `is`.
        """

        if isinstance(other, dict):
            other = Molecule(**other)
        elif isinstance(other, Molecule):
            pass
        else:
            raise TypeError(
                f"Comparison molecule not understood of type '{type(other)}'."
            )

        return self.get_hash() == other.get_hash()

    def pretty_print(self):
        """Prints the molecule. Not sure yet what this is used for, but I have modified the
        original implementation in qcelemental.
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
            for j in range(3):
                text += """  {0:17.12f}""".format(self.geometry[i][j])
            text += "\n"

        return text

    @property
    def hash_fields(self):
        return [
            "symbols",
            "masses",
            "molecular_charge",
            "molecular_multiplicity",
            "real",
            "geometry",
            "velocities",
            "forces",
            "fragments",
            "fragment_charges",
            "fragment_multiplicities",
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
                if field == "geometry":
                    data = qcelemental.models.molecule.float_prep(data, GEOMETRY_NOISE)
                elif field == "velocities":
                    data = qcelemental.models.molecule.float_prep(data, VELOCITY_NOISE)
                elif field == "forces":
                    data = qcelemental.models.molecule.float_prep(data, FORCE_NOISE)
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
        **kwargs,
    ) -> "Molecule":
        """
        Constructs a Molecule object from a file.
        Parameters
        ----------
        filename : str
            The atomic positions filename to read
        top_filename: str, optional
            The topology i.e. connectivity filename to read
        dtype : str, optional
            The type of file to interpret. If not set, mmelemental attempts to discover the file type.
        **kwargs
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

            import json

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

        tkmol = cls._tkFromFile(
            filename=filename, top_filename=top_filename, dtype=dtype
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
        **kwargs
            Additional kwargs to pass to the constructors.
        Returns
        -------
        Molecule
            A constructed Molecule class.
        """
        if isinstance(data, str):
            if not dtype:
                raise ValueError(
                    "You must supply dtype for proper interpretation of symbolic data. See the :class:``Identifiers`` class."
                )
            try:
                code = ChemCode(code=data, dtype=dtype)
                symbols = list(code.code)  # this is garbage, must be replaced
                return Molecule(identifiers={dtype: data}, symbols=symbols)
            except Exception as e:
                raise ValueError(
                    f"Failed in interpreting {data} as a valid code. Exception: {e}"
                )
        elif isinstance(data, dict):
            kwargs.update(data)
            return cls(**kwargs)

        return data.to_schema(**kwargs)

    def to_file(
        self,
        filename: str,
        dtype: Optional[str] = None,
        translator: Optional[str] = None,
        **kwargs,
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
            Translator name e.g. mmic_parmed. If unset, MMElemental attempts to find an appropriate 
            translator if it is registered in the :class:``TransComponent`` class. 
        **kwargs
            Additional kwargs to pass to the constructor.
        """
        if not dtype:
            from pathlib import Path

            ext = Path(filename).suffix
        else:
            ext = "." + dtype

        qcelemental_handle = ext in qcelemental.models.molecule._extension_map

        if qcelemental_handle:
            dtype = qcelemental.models.molecule._extension_map[ext]
            super().to_file(
                filename, dtype, **kwargs
            )  # qcelemental handles serialization
        else:  # look for an installed mmic_translator
            translator = TransComponent.find_molwrite_tk(ext)
            tkmol = self._tkFromSchema(translator=translator, **kwargs)
            tkmol.to_file(filename, dtype=dtype, **kwargs)  # pass dtype?

    def to_data(
        self,
        dtype: str,
        translator: Optional[str] = None,
        **kwargs: Optional[Dict[str, Any]],
    ) -> ToolkitModel:
        """Converts Molecule to toolkit-specific molecule (e.g. rdkit, MDAnalysis, parmed).
        Parameters
        ----------
        dtype: str
            The type of data object to convert to e.g. MDAnalysis, rdkit, parmed, etc.
        translator: Optional[str], optional
            Translator name e.g. mmic_rdkit. If unset, MMElemental attempts to find an appropriate 
            translator if it is registered in the :class:``TransComponent`` class. 
        **kwargs
            Additional kwargs to pass to the constructor.
        Returns
        -------
        ToolkitModel
            Toolkit-specific molecule model
        """
        return self._tkFromSchema(dtype=dtype, translator=translator, **kwargs)

    def _tkFromSchema(
        self,
        dtype: str = None,
        translator: Optional[str] = None,
        **kwargs: Optional[Dict[str, Any]],
    ) -> ToolkitModel:
        """Helper function that constructs a toolkit-specific molecule from MMSchema molecule.
        Which toolkit-specific component is called depends on which package is installed on the system.
        Parameters
        ----------
        dtype: Optional[str], optional
            The type of data object to convert to e.g. MDAnalysis, rdkit, parmed, etc.
        translator: Optional[str], optional
            Translator name e.g. mmic_parmed. If unset, MMElemental attempts to find an appropriate 
            translator if it is registered in the :class:``TransComponent`` class. 
        **kwargs
            Additional kwargs to pass to the constructor.
        Returns
        -------
        ToolkitModel
            Toolkit-specific molecule model
        """

        if not translator:
            if not dtype:
                raise ValueError(
                    f"Either translator or dtype must be supplied when calling {__name__}."
                )
            translator = TransComponent.find_trans(dtype)

        if translator == "qcelemental":
            qmol = qcelemental.models.molecule.Molecule.to_data(
                data, orient=False, validate=False, **kwargs
            )
            return Molecule(orient=orient, validate=validate, **qmol.to_dict())
        elif importlib.util.find_spec(translator):
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

    @classmethod
    def _tkFromFile(
        cls,
        filename: str,
        top_filename: Optional[str] = None,
        dtype: Optional[str] = None,
        translator: Optional[str] = None,
        **kwargs: Optional[Dict[str, Any]],
    ) -> ToolkitModel:
        """Helper function that constructs a toolkit-specific molecule from an input file.
        Which toolkit-specific component is called depends on which package is installed on the system.
        Parameters
        ----------
        filename: str
            Molecule file name. The file should define the molecule.
        top_filename: Optional[str], optional
            Topology file name. The file should define the molecular connectivity.
        dtype: Optional[str], optional
            The type of data file object to read e.g. gro, pdb, etc.
        translator: Optional[str], optional
            Translator name e.g. mmic_mda. If unset, MMElemental attempts to find an appropriate 
            translator if it is registered in the :class:``TransComponent`` class. 
            Additional kwargs to pass to the constructor.
        Returns
        -------
        ToolkitModel
        """

        fileobj = FileOutput(path=filename) if filename else None
        top_fileobj = FileOutput(path=top_filename) if top_filename else None

        dtype = dtype or fileobj.ext

        if not translator:
            from mmelemental.components.trans.supported import reg_trans

            reg_trans = list(reg_trans)

            while not translator:
                translator = TransComponent.find_molread_tk(dtype, trans=reg_trans)
                if not translator:
                    raise ValueError(
                        f"Could not read xyz file with ext {dtype}. Please install an appropriate translator."
                    )
                # Make sure we can import the translator module
                if importlib.util.find_spec(translator):
                    mod = importlib.import_module(translator)

                # If top if supplied, make sure the translator supports the top file extension
                if top_fileobj:
                    top_dtype = top_fileobj.ext
                    if top_dtype not in mod.ffread_ext_maps:
                        reg_trans.remove(translator)
                        translator = None
                    if not len(reg_trans):
                        raise ValueError(
                            f"Could not read xyz and top files with exts {dtype} and {top_dtype}. \
                            Please install an appropriate translator."
                        )

        tkmol = mod._classes_map.get("Molecule")

        if not tkmol:
            raise ValueError(
                f"No Molecule model found while looking in translator: {translator}."
            )

        return tkmol.from_file(
            filename=fileobj.abs_path if fileobj else None,
            top_filename=top_fileobj.abs_path if top_fileobj else None,
            dtype=dtype.strip("."),
        )
