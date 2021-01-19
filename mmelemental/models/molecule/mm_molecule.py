import qcelemental
from qcelemental.models.types import Array
from typing import List, Tuple, Optional, Any, Dict, Union
from pydantic import validator, Field, ValidationError
from mmelemental.components.io.molreader_component import TkMolReaderComponent
from mmelemental.models.molecule.io_molecule import MolInput, MolOutput
from mmelemental.models.chem.codes import ChemCode
from mmelemental.models.util.input import FileInput
from mmelemental.models.util.output import FileOutput
from mmic.components.blueprints.generic_component import GenericComponent
from mmelemental.models.base import Nothing

class Identifiers(qcelemental.models.molecule.Identifiers):
    """ 
    An extension of the qcelemental.models.molecule.Identifiers for RDKit constructors.
    See `link <https://rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html>`_ for more info. 
    """
    smiles: Optional[ChemCode] = Field(
        None,
        description="A simplified molecular-input line-entry system code."
    )
    smarts: Optional[ChemCode] = Field(
        None,
        description="A SMILES arbitrary target specification code for defining substructures."
    )
    inchi: Optional[ChemCode] = Field(
        None,
        description="An international chemical identifier code."
    )
    sequence: Optional[ChemCode] = Field(
        None,
        description="A sequence code from RDKit (currently only supports peptides)."
    )
    fasta: Optional[ChemCode] = Field(
        None,
        description="A FASTA code (currently only supports peptides)."
    )
    helm: Optional[ChemCode] = Field(
        None,
        description="A HELM code (currently only supports peptides)."
    )

class Molecule(qcelemental.models.Molecule):
    """
    An MMSchema representation of a Molecule based on QCSchema. This model contains data for symbols, geometry, 
    connectivity, charges, residues, etc. while also supporting a wide array of I/O and manipulation capabilities.
    Molecule objects geometry, masses, and charges are truncated to 8, 6, and 4 decimal places respectively 
    to assist with duplicate detection.
    """
    symbols: Array[str] = Field(
        None,
        description = "An ordered (natom,) array-like object of atomic elemental symbols. The index of "
        "this attribute sets atomic order for all other per-atom setting like ``real`` and the first "
        "dimension of ``geometry``. Ghost/Virtual atoms must have an entry in this array-like and are "
        "indicated by the matching the 0-indexed indices in ``real`` field.",
    )
    geometry: Optional[Array[float]] = Field(
        None,
        description="An ordered (natom,3) array-like for XYZ atomic positions in Angstrom. \n
        "Can also accept arrays which can be mapped to (natom,3) such as a 1-D list of length 3*natom, "
        "or the serialized version of the array in (3*natom,) shape; all forms will be reshaped to "
        "(natom,3) for this attribute.",
    )
    velocities: Optional[Array[float]] = Field(
        None,
        description="An ordered (natoms,3) array-like for XYZ atomic velocities in Angstrom/ps. \n
        "Can also accept arrays which can be mapped to (natoms,3) such as a 1-D list of length 3*natoms, "
        "or the serialized version of the array in (3*natoms,) shape; all forms will be reshaped to "
        "(natoms,3) for this attribute.",
    )
    forces: Optional[Array[float]] = Field(
        None,
        description="An ordered (natoms,3) array-like for XYZ atomic velocities in kJ/mol*Angstrom. \n
        "Can also accept arrays which can be mapped to (natoms,3) such as a 1-D list of length 3*natoms, "
        "or the serialized version of the array in (3*natoms,) shape; all forms will be reshaped to "
        "(natoms,3) for this attribute.",
    )
    angles: Optional[List[Tuple[int, int, int]]] = Field(
        None,
        description = "Bond angles in degrees for three connected atoms."
    )
    dihedrals: Optional[List[Tuple[int, int, int, int, int]]] = Field(
        None,
        description = 'Dihedral/torsion angles in degrees between planes through two sets of three atoms, having two atoms in common.')
    improper_dihedrals: Optional[List[Tuple[int, int, int, int, int]]] = Field(
        None,
        description = 'Improper dihedral/torsion angles in degrees between planes through two sets of three atoms, having two atoms in common.')
    residues: Optional[List[Tuple[str, int]]] = Field(
        None, 
        description = "A list of (residue_name, residue_num) of connected atoms constituting the building block (monomer) "
        "of a polymer. Order follows atomic indices from 0 till Natoms-1. Residue number starts from 1."
        "\n"
        "E.g. ('ALA', 1) means atom 0 belongs to aminoacid alanine with residue number 1."
        )
    chains: Optional[Dict[str, List[int]]] = Field(
        None, description = "A sequence of connected residues (i.e. polymers) forming a subunit that is not bonded to any "
        "other subunit. For example, a hemoglobin molecule consists of four chains that are not connected to one another."
    )
    segments: Optional[Dict[str, List[int]]] = Field(
        None, 
        description = "..."
    )
    names: Optional[List[str]] = Field(
        None, 
        description = "A list of atomic label names."
    )
    identifiers: Optional[Identifiers] = Field(
        None, 
        description = "An optional dictionary of additional identifiers by which this Molecule can be referenced, "
        "such as INCHI, SMILES, SMARTS, etc. See the :class:``Identifiers`` model for more details."
    )
    rotateBonds: Optional[List[Tuple[int, int]]] = Field(
        None, 
        description = "A list of bonded atomic indices: (atom1, atom2), specifying rotatable bonds in the molecule."
    )
    rigidBonds: Optional[List[Tuple[int, int]]] = Field(
        None, description = "A list of bonded atomic indices: (atom1, atom2), specifying rigid bonds in the molecule."
    )

    # Constructors
    @classmethod
    def from_file(cls, filename: Union[FileInput, str], top: Union[FileInput, str] = None, dtype: Optional[str] = None, 
        *, orient: bool = False, **kwargs) -> "Molecule":
        """
        Constructs a Molecule object from a file.
        Parameters
        ----------
        filename : str
            The coords filename to read
        top: str
            The topology i.e. connectivity filename to read
        dtype : str, optional
            The type of file to interpret. If not set, mmelemental attempts to discover the file type.
        orient : bool, optional
            Orientates the molecule to a standard frame or not.
        **kwargs
            Any additional keywords to pass to the constructor
        Returns
        -------
        Molecule
            A constructed Molecule class.
        """
        if not isinstance(filename, FileInput):
            fileobj = FileInput(path=filename, dtype=dtype)
        else:
            fileobj = FileInput(path=filename.path, dtype=dtype)

        if top and not isinstance(top, FileInput):
            top = FileInput(path=top, dtype=dtype)
        elif top and isinstance(top, FileInput):
            top = FileInput(path=top.path, dtype=dtype)
 
        if top:
            mol_input = MolInput(file=fileobj, top_file=top)
        else:
            mol_input = MolInput(file=fileobj)

        mol = TkMolReaderComponent.compute(mol_input)

        return cls.from_data(mol, dtype=mol.dtype)
        
    @classmethod
    def from_data(cls, data: Any, dtype: Optional[str] = None, *,
        orient: bool = False, validate: bool = None, **kwargs: Dict[str, Any]) -> "Molecule":
        """
        Constructs a Molecule object from a data object.
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
        Molecule
            A constructed Molecule class.
        """
        if isinstance(data, str):
            try:
                code = ChemCode(code=data)
                mol_input = MolInput(code=code, args={'validate': validate, 'orient': orient, 'kwargs': kwargs})
            except:
                raise ValueError
        elif isinstance(data, ChemCode):
            mol_input = MolInput(code=data, args={'validate': validate, 'orient': orient, 'kwargs': kwargs})
        else:
            # Let's hope this is a toolkit-specific molecule and pass it as data
            mol_input = MolInput(data=data, args={'validate': validate, 'orient': orient, 'kwargs': kwargs})
        
        return MolReaderComponent.compute(mol_input)

    def to_file(self, filename: str, dtype: Optional[str] = None, mode: str = 'w') -> Nothing:
        """ Writes the Molecule to a file.
        Parameters
        ----------
        filename : str
            The filename to write to
        dtype : Optional[str], optional
            The type of file to write, attempts to infer dtype from the filename if not provided.
        mode: str
            Write new file or overwrite existing file (w) or append (a) to existing file.
        """
        fileobj = FileOutput(path=filename, dtype=dtype)
        mol_input = MolOutput(file=fileobj, mol=self, mode=mode)
        toolkit = mol_input.files_toolkit()

        if toolkit == 'qcelemental':
            super().to_file(filename, dtype, mode)
        else:
            MolWriterComponent.compute(mol_input)

    def to_data(self, dtype: str) -> "ToolkitMolecule":
        """ Converts Molecule to toolkit-specific molecule (e.g. rdkit). """

        if dtype == 'rdkit':
            from mmelemental.components.trans.rdkit_component import MolToRDKitComponent
            return MolToRDKitComponent.compute(self).mol
        else:
            raise NotImplementedError(f'Data type {dtype} not available.')


class MolReaderComponent(GenericComponent):
    """ Factory component that constructs a Molecule object from MolInput.
    Which toolkit-specific component is called depends on data type and 
    which toolkits are installed on the system."""

    @classmethod
    def input(cls):
        return MolInput

    @classmethod
    def output(cls):
        return Molecule

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None) -> Tuple[bool, Dict[str, Any]]:

        if inputs.args:
            orient = inputs.args.get('orient')
            validate = inputs.args.get('validate')
            kwargs = inputs.args.get('kwargs')
        else:
            orient, validate, kwargs = False, None, None

        if inputs.data:
            dtype = inputs.data.dtype
            if dtype == 'qcelemental':
                qmol = qcelemental.models.molecule.Molecule.from_data(data, dtype, orient=orient, validate=validate, **kwargs)
                return True, Molecule(orient=orient, validate=validate, **qmol.to_dict())
            elif dtype == 'rdkit':
                from mmelemental.components.trans.rdkit_component import RDKitToMolComponent
                return True, RDKitToMolComponent.compute(inputs)
            elif dtype == 'parmed':
                from mmelemental.components.trans.parmed_component import ParmedToMolComponent
                return True, ParmedToMolComponent.compute(inputs)               
            else:
                raise NotImplementedError(f'Data type not yet supported: {dtype}.')
        # Only RDKit is handling chem codes and file objects for now!
        elif inputs.code:
            from mmelemental.components.trans.rdkit_component import RDKitToMolComponent
            return True, RDKitToMolComponent.compute(inputs)
        elif inputs.file:
            toolkit = inputs.files_toolkit()
            if toolkit == 'rdkit':
                from mmelemental.components.trans.rdkit_component import RDKitToMolComponent
                return True, RDKitToMolComponent.compute(inputs)
            elif toolkit == 'parmed':
                from mmelemental.components.trans.parmed_component import ParmedToMolComponent
                return True, ParmedToMolComponent.compute(inputs)
        else:
            raise NotImplementedError('Molecules can be instantiated from codes, files, or other data objects.')

class MolWriterComponent(GenericComponent):
    """ Factory component that constructs a Molecule object from MolInput.
    Which toolkit-specific component is called depends on MolInput.data.dtype."""

    @classmethod
    def input(cls):
        return MolOutput

    @classmethod
    def output(cls):
        return Nothing

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None) -> Tuple[bool, None]:
        
        if isinstance(inputs, dict):
            inputs = MolReaderComponent.input()(**inputs)

        if inputs.args:
            orient = inputs.args.get('orient')
            validate = inputs.args.get('validate')
            kwargs = inputs.args.get('kwargs')
        else:
            orient, validate, kwargs = False, None, None

        toolkit = inputs.files_toolkit()
        filename = inputs.file.abs_path
        dtype = inputs.file.dtype or inputs.file.ext
        mode = inputs.file.mode

        if toolkit == 'rdkit':
            from mmelemental.components.trans.rdkit_component import MolToRDKitComponent
            from rdkit import Chem
            
            rdkmol = MolToRDKitComponent.compute(inputs.mol)
            if mode != 'w':
                raise NotImplementedError('rdkit supports only write mode for file output. Ouch!')

            if dtype == '.pdb':
                writer = Chem.PDBWriter(filename)
            elif dtype == '.sdf':
                writer = Chem.SDWriter(fp)
            elif dtype == '.smi':
                writer = Chem.SmilesWriter(fp)
            else:
                raise NotImplementedError(f'File format {dtype} not supported by rdkit.')

            writer.write(rdkmol.mol)
            writer.close()

        elif toolkit == 'parmed':
            overwrite = True if inputs.file.mode == 'w' else False
            from mmelemental.components.trans.parmed_component import MolToParmedComponent
            pmol = MolToParmedComponent.compute(inputs.mol)
            #print([atom.name for atom in pmol.mol.atoms])
            #print([(residue.name, residue.atoms) for residue in pmol.mol.residues])
            pmol.mol.save(filename, overwrite=overwrite)
        else:
            raise ValueError(f'Data type not yet supported: {dtype}')

        return True, Nothing()