import qcelemental
from qcelemental.models.types import Array
from typing import List, Tuple, Optional, Any, Dict, Union
from pydantic import Field, constr
import importlib

# Import MM models
from mmelemental.models.molecule.io_mol import MolInput, MolOutput
from mmelemental.models.util.output import FileOutput
from mmelemental.models.base import ToolkitModel
from mmelemental.models.chem.codes import ChemCode
from mmelemental.models.base import Provenance, provenance_stamp

# Import MM components
from mmelemental.components.trans import TransComponent
from mmic.components.blueprints.generic_component import GenericComponent


__all__ = ["Molecule"]


class Identifiers(qcelemental.models.molecule.Identifiers):
    """
    An extension of the qcelemental.models.molecule.Identifiers for RDKit constructors.
    See `link <https://rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html>`_ for more info.
    """

    smiles: Optional[ChemCode] = Field(
        None, description="A simplified molecular-input line-entry system code."
    )
    smarts: Optional[ChemCode] = Field(
        None,
        description="A SMILES arbitrary target specification code for defining substructures.",
    )
    inchi: Optional[ChemCode] = Field(
        None, description="An international chemical identifier code."
    )
    sequence: Optional[ChemCode] = Field(
        None,
        description="A sequence code from RDKit (currently only supports peptides).",
    )
    fasta: Optional[ChemCode] = Field(
        None, description="A FASTA code (currently only supports peptides)."
    )
    helm: Optional[ChemCode] = Field(
        None, description="A HELM code (currently only supports peptides)."
    )


mmschema_molecule_default = "mmschema_molecule"


class Molecule(qcelemental.models.Molecule):
    """
    A representation of a Molecule in MM based on QCSchema. This model contains data for symbols, geometry,
    connectivity, charges, residues, etc. while also supporting a wide array of I/O and manipulation capabilities.
    Molecule objects geometry, masses, and charges are truncated to 8, 6, and 4 decimal places respectively
    to assist with duplicate detection.
    """

    schema_name: constr(strip_whitespace=True, regex="^(mmschema_molecule)$") = Field(  # type: ignore
        mmschema_molecule_default,
        description=(
            f"The MMSchema specification to which this model conforms. Explicitly fixed as {mmschema_molecule_default}."
        ),
    )
    symbols: Optional[Array[str]] = Field(
        None,
        description="An ordered (natom,) array-like object of atomic elemental symbols. The index of "
        "this attribute sets atomic order for all other per-atom setting like ``real`` and the first "
        "dimension of ``geometry``. Ghost/Virtual atoms must have an entry in this array-like and are "
        "indicated by the matching the 0-indexed indices in ``real`` field.",
    )
    masses_units: Optional[str] = Field(
        "amu",
        description="Units for atomic masses. Defaults to unified atomic mass unit.",
    )
    molecular_charge_units: Optional[str] = Field(
        "eV", description="Units for molecular charge. Defaults to electron Volt."
    )
    geometry: Optional[Array[float]] = Field(
        None,
        description="An ordered (natom,3) array-like for XYZ atomic positions in Angstrom. "
        "Can also accept arrays which can be mapped to (natom,3) such as a 1-D list of length 3*natom, "
        "or the serialized version of the array in (3*natom,) shape; all forms will be reshaped to "
        "(natom,3) for this attribute. Default unit is Angstroms.",
    )
    geometry_units: Optional[str] = Field(
        "angstrom", description="Units for atomic geometry. Defaults to Angstroms."
    )
    velocities: Optional[Array[float]] = Field(
        None,
        description="An ordered (natoms,3) array-like for XYZ atomic velocities in Angstrom/ps. "
        "Can also accept arrays which can be mapped to (natoms,3) such as a 1-D list of length 3*natoms, "
        "or the serialized version of the array in (3*natoms,) shape; all forms will be reshaped to "
        "(natoms,3) for this attribute. Default unit is Angstroms/femtoseconds.",
    )
    velocities_units: Optional[str] = Field(
        "angstrom/fs",
        description="Units for atomic velocities. Defaults to Angstroms/femtoseconds.",
    )
    forces: Optional[Array[float]] = Field(
        None,
        description="An ordered (natoms,3) array-like for XYZ atomic velocities in kJ/mol*Angstrom. "
        "Can also accept arrays which can be mapped to (natoms,3) such as a 1-D list of length 3*natoms, "
        "or the serialized version of the array in (3*natoms,) shape; all forms will be reshaped to "
        "(natoms,3) for this attribute. Default unit is KiloJoules/mol.Angstroms.",
    )
    forces_units: Optional[str] = Field(
        "kJ/(mol*angstrom)",
        description="Units for atomic forces. Defaults to KiloJoules/mol.Angstroms",
    )
    angles: Optional[List[Tuple[int, int, int]]] = Field(
        None, description="Bond angles in degrees for three connected atoms."
    )
    dihedrals: Optional[List[Tuple[int, int, int, int, int]]] = Field(
        None,
        description="Dihedral/torsion angles in degrees between planes through two sets of three atoms, having two atoms in common.",
    )
    improper_dihedrals: Optional[List[Tuple[int, int, int, int, int]]] = Field(
        None,
        description="Improper dihedral/torsion angles in degrees between planes through two sets of three atoms, having two atoms in common.",
    )
    residues: Optional[List[Tuple[str, int]]] = Field(
        None,
        description="A list of (residue_name, residue_num) of connected atoms constituting the building block (monomer) "
        "of a polymer. Order follows atomic indices from 0 till Natoms-1. Residue number starts from 1."
        "\n"
        "E.g. ('ALA', 1) means atom 0 belongs to aminoacid alanine with residue number 1.",
    )
    chains: Optional[Dict[str, List[int]]] = Field(
        None,
        description="A sequence of connected residues (i.e. polymers) forming a subunit that is not bonded to any "
        "other subunit. For example, a hemoglobin molecule consists of four chains that are not connected to one another.",
    )
    segments: Optional[Dict[str, List[int]]] = Field(None, description="...")
    names: Optional[Union[List[str], Array[str]]] = Field(
        None, description="A list of atomic label names."
    )
    identifiers: Optional[Identifiers] = Field(
        None,
        description="An optional dictionary of additional identifiers by which this Molecule can be referenced, "
        "such as INCHI, SMILES, SMARTS, etc. See the :class:``Identifiers`` model for more details.",
    )
    rotateBonds: Optional[List[Tuple[int, int]]] = Field(
        None,
        description="A list of bonded atomic indices: (atom1, atom2), specifying rotatable bonds in the molecule.",
    )
    rigidBonds: Optional[List[Tuple[int, int]]] = Field(
        None,
        description="A list of bonded atomic indices: (atom1, atom2), specifying rigid bonds in the molecule.",
    )
    provenance: Provenance = Field(
        provenance_stamp(__name__),
        description="The provenance information about how this object (and its attributes) were generated, "
        "provided, and manipulated.",
    )

    def __init__(self, **kwargs: Any) -> None:
        """Initializes the molecule object from dictionary-like values.
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

    # Constructors
    @classmethod
    def from_file(
        cls,
        filename: Optional[str] = None,
        top_filename: Optional[str] = None,
        dtype: Optional[str] = None,
        **kwargs,
    ) -> "Molecule":
        """
        Constructs a Molecule object from a file.
        Parameters
        ----------
        filename : str, optional
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
        if not top_filename and not filename:
            raise TypeError(
                "You must supply at least one of the following: filename or top_filename."
            )

        tkmol = cls._tkFromFile(
            filename=filename, top_filename=top_filename, dtype=dtype
        )
        return cls.from_data(tkmol, dtype=tkmol.dtype)

    @classmethod
    def from_data(
        cls,
        data: Any,
        dtype: Optional[str] = None,
        **kwargs: Dict[str, Any],
    ) -> "Molecule":
        """
        Constructs a Molecule object from a data object.
        Parameters
        ----------
        data: Any
            Data to construct Molecule from
        dtype: str, optional
            How to interpret the data, if not passed attempts to discover this based on input type.
        **kwargs
            Additional kwargs to pass to the constructors. kwargs take precedence over data.
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
                data = ChemCode(code=data)
                return Molecule(identifiers={dtype: data})
            except Exception:
                raise ValueError("Failed in interpreting {data} as {}.")

        return data.to_schema(**kwargs)

    def to_file(self, filename: str, dtype: Optional[str] = None, **kwargs) -> None:
        """Writes the Molecule to a file.
        Parameters
        ----------
        filename : str
            The filename to write to
        dtype : Optional[str], optional
            The type of file to write (e.g. pdb, gro, etc.), attempts to infer dtype from
            file extension if not provided.
        **kwargs
            Additional kwargs to pass to the constructor.
        """
        if not dtype:
            from pathlib import Path

            ext = Path(filename).suffix

        qcelemental_handle = ext in qcelemental.models.molecule._extension_map

        if qcelemental_handle:
            super().to_file(
                filename, ext, **kwargs
            )  # replace this with openbabel that can handle xyz files
        else:  # look for an installed mmic_translator
            translator = TransComponent.find_molwrite_tk(ext)
            tkmol = self._tkFromSchema(translator=translator, **kwargs)
            tkmol.to_file(filename, dtype=ext.strip("."), **kwargs)  # pass dtype?

    def to_data(self, dtype: str, **kwargs) -> ToolkitModel:
        """Converts Molecule to toolkit-specific molecule (e.g. rdkit, MDAnalysis, parmed).
        Parameters
        ----------
        dtype: str
            The type of data object to convert to e.g. MDAnalysis, rdkit, parmed, etc.
        **kwargs
            Additional kwargs to pass to the constructor.
        Returns
        -------
        ToolkitModel
            Toolkit-specific molecule model
        """
        return self._tkFromSchema(dtype=dtype, **kwargs)

    def _tkFromSchema(
        self, translator: str = None, dtype: str = None, **kwargs
    ) -> ToolkitModel:
        """Helper function that constructs a toolkit-specific molecule from MMSchema molecule.
        Which toolkit-specific component is called depends on which package is installed on the system.
        Parameters
        ----------
        translator: str
            translator name e.g. mmic_parmed. The translator should be registered in the :class:``TransComponent`` class.
        dtype: str
            The type of data object to convert to e.g. MDAnalysis, rdkit, parmed, etc.
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
        cls, filename: str = None, top_filename: str = None, dtype: str = None, **kwargs
    ) -> ToolkitModel:
        """Helper function that constructs a toolkit-specific molecule from an input file.
        Which toolkit-specific component is called depends on which package is installed on the system.
        Parameters
        ----------
        filename: str, optional
            Molecule file name. The file should define the molecule.
        top_filename: str, optional
            Topology file name. The file should define the molecular connectivity.
        dtype: str
            The type of data file object to read e.g. gro, pdb, etc.
        **kwargs
            Additional kwargs to pass to the constructor.
        Returns
        -------
        ToolkitModel
        """

        fileobj = FileOutput(path=filename) if filename else None
        top_fileobj = FileOutput(path=top_filename) if top_filename else None

        if top_filename:
            dtype = dtype or top_fileobj.ext
        elif filename:
            dtype = dtype or fileobj.ext
        else:
            raise TypeError(
                "At least one file must be supplied! You must define either the connectivity or molecule, or both."
            )
        translator = TransComponent.find_molread_tk(dtype)

        if not translator:
            raise ValueError(
                f"Could not read file with ext {dtype}. Please install an appropriate translator."
            )

        if importlib.util.find_spec(translator):
            mod = importlib.import_module(translator)
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
