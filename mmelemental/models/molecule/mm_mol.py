import qcelemental
from qcelemental.models.types import Array
from typing import List, Tuple, Optional, Any, Dict, Union
from pydantic import validator, Field, ValidationError
from mmelemental.components.io.molreader_component import TkMolReaderComponent
from mmelemental.models.molecule.io_mol import MolInput, MolOutput, Translators
from .gen_mol import ToolkitMol
from mmelemental.models.chem.codes import ChemCode
from mmelemental.models.util.input import FileInput
from mmelemental.models.util.output import FileOutput
from mmic.components.blueprints.generic_component import GenericComponent
from mmelemental.models.base import Nothing
import importlib


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


class Mol(qcelemental.models.Molecule):
    """
    A representation of a Molecule in MM based on QCSchema. This model contains data for symbols, geometry,
    connectivity, charges, residues, etc. while also supporting a wide array of I/O and manipulation capabilities.
    Molecule objects geometry, masses, and charges are truncated to 8, 6, and 4 decimal places respectively
    to assist with duplicate detection.
    """

    symbols: Array[str] = Field(
        None,
        description="An ordered (natom,) array-like object of atomic elemental symbols. The index of "
        "this attribute sets atomic order for all other per-atom setting like ``real`` and the first "
        "dimension of ``geometry``. Ghost/Virtual atoms must have an entry in this array-like and are "
        "indicated by the matching the 0-indexed indices in ``real`` field.",
    )
    geometry: Optional[Array[float]] = Field(
        None,
        description="An ordered (natom,3) array-like for XYZ atomic positions in Angstrom. "
        "Can also accept arrays which can be mapped to (natom,3) such as a 1-D list of length 3*natom, "
        "or the serialized version of the array in (3*natom,) shape; all forms will be reshaped to "
        "(natom,3) for this attribute.",
    )
    velocities: Optional[Array[float]] = Field(
        None,
        description="An ordered (natoms,3) array-like for XYZ atomic velocities in Angstrom/ps. "
        "Can also accept arrays which can be mapped to (natoms,3) such as a 1-D list of length 3*natoms, "
        "or the serialized version of the array in (3*natoms,) shape; all forms will be reshaped to "
        "(natoms,3) for this attribute.",
    )
    forces: Optional[Array[float]] = Field(
        None,
        description="An ordered (natoms,3) array-like for XYZ atomic velocities in kJ/mol*Angstrom. "
        "Can also accept arrays which can be mapped to (natoms,3) such as a 1-D list of length 3*natoms, "
        "or the serialized version of the array in (3*natoms,) shape; all forms will be reshaped to "
        "(natoms,3) for this attribute.",
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

    # Constructors
    @classmethod
    def from_file(
        cls,
        filename: Optional[str] = None,
        top_filename: Optional[str] = None,
        dtype: Optional[str] = None,
        *,
        orient: bool = False,
        **kwargs,
    ) -> "Mol":
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
        orient : bool, optional
            Orientates the molecule to a standard frame or not.
        **kwargs
            Any additional keywords to pass to the constructor
        Returns
        -------
        Mol
            A constructed Mol class.
        """
        if top_filename and filename:
            mol_input = MolInput(file=filename, top_file=top_filename, dtype=dtype)
        elif filename:
            mol_input = MolInput(file=filename, dtype=dtype)
        elif top_filename:
            mol_input = MolInput(file=filename, dtype=dtype)
        else:
            raise TypeError(
                "You must supply at least one of the following: filename or top_filename."
            )

        tkmol = TkMolReaderComponent.compute(mol_input)
        return cls.from_data(tkmol, dtype=tkmol.dtype)

    @classmethod
    def from_data(
        cls,
        data: Any,
        dtype: Optional[str] = None,
        *,
        orient: bool = False,
        validate: bool = None,
        **kwargs: Dict[str, Any],
    ) -> "Mol":
        """
        Constructs a Mol object from a data object.
        Parameters
        ----------
        data: Any
            Data to construct Molecule from
        dtype: str, optional
            How to interpret the data, if not passed attempts to discover this based on input type.
        orient: bool, optional
            Orientates the molecule to a standard frame or not.
        validate: bool, optional
            Validates the molecule or not.
        **kwargs
            Additional kwargs to pass to the constructors. kwargs take precedence over data.
        Returns
        -------
        Mol
            A constructed Mol class.
        """
        if isinstance(data, str):
            try:
                data = ChemCode(code=data)
            except:
                raise ValueError

        return data.to_data(orient=orient, validate=validate, **kwargs)

    def to_file(
        self, filename: str, dtype: Optional[str] = None, mode: str = "w", **kwargs
    ) -> Nothing:
        """Writes the Molecule to a file.
        Parameters
        ----------
        filename : str
            The filename to write to
        dtype : Optional[str], optional
            The type of file to write, attempts to infer dtype from the filename if not provided.
        mode: str
            Write new file or overwrite existing file (w) or append (a) to existing file.
        **kwargs
            Additional kwargs to pass to the constructor. kwargs take precedence over data.
        """
        inputs = MolOutput(file=filename, mol=self, mode=mode)
        # tkmol = TkMolWriterComponent.compute(inputs)
        # tkmol.to_file(filename, dtype, mode, **kwargs)

    def to_data(self, dtype: str, **kwargs) -> ToolkitMol:
        """Converts Molecule to toolkit-specific molecule (e.g. rdkit, MDAnalysis, parmed).
        Parameters
        ----------
        dtype: str
            The type of data object to convert to.
        **kwargs
            Additional kwargs to pass to the constructor.
        """
        inputs = MolOutput(mol=self, dtype=dtype, kwargs=kwargs)
        return FromMolComponent.compute(inputs)


class FromMolComponent(GenericComponent):
    """Factory component that reads a Mol object and constructs a toolkit-specific molecule.
    Which toolkit-specific component is called depends on which package is installed on the system."""

    @classmethod
    def input(cls):
        return MolOutput

    @classmethod
    def output(cls):
        return ToolkitMol

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, None]:

        import inspect

        translator = Translators.find_trans(inputs.dtype)

        if translator == "mmic_qcelemental":
            qmol = qcelemental.models.molecule.Molecule.to_data(
                data, orient=orient, validate=validate, **inputs.kwargs
            )
            return True, Mol(orient=orient, validate=validate, **qmol.to_dict())
        elif importlib.util.find_spec(translator):
            mod = importlib.import_module(translator + ".models")
            models = inspect.getmembers(mod, inspect.isclass)

            tkmol = [mol for _, mol in models if issubclass(mol, ToolkitMol)]

            if len(tkmol) > 1:
                raise ValueError(
                    "More than 1 compatible model found:"
                    + (" {} " * len(tkmol)).format(*tkmol)
                )
            elif not tkmol:
                raise ValueError(
                    f"No compatible model found while looking in translator: {translator}."
                )

            tkmol = tkmol[0]
            return True, tkmol.from_data(inputs.mol)
        else:
            raise NotImplementedError(
                f"Translator for {inputs.dtype} not yet available."
            )
