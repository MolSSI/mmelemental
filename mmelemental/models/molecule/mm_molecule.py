import qcelemental
from qcelemental.models.types import Array
from typing import List, Tuple, Optional, Any, Dict, Union
import os, sys
import random
import string
import numpy
from pydantic import validator, Field, ValidationError
from mmelemental.components.molreader_component import TkMoleculeReaderComponent, MoleculeReaderComponent
from mmelemental.models.molecule.mol_reader import MoleculeReaderInput
from mmelemental.models.chem.codes import ChemCode
from mmelemental.models.util.input import FileInput
from pathlib import Path

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
    coordinates: Array[float] = Field(
        None,
        description = "An ordered (natoms, 3) array-like for XYZ atomic coordinates."
    )
    velocities: Array[float] = Field(
        None,
        description = "An ordered (natoms, 3) array-like for XYZ atomic velocities",
    )
    forces: Array[float] = Field(
        None,
        description = "An ordered (natoms, 3) array-like for XYZ atomic forces",
    )
    angles: Optional[List[Tuple[int, int, int]]] = Field(
        None,
        description = "Bond angles of three connected atoms."
    )
    dihedrals: Optional[List[Tuple[int, int, int, int, int]]] = Field(
        None,
        description = 'Dihedral/torsion angles between planes through two sets of three atoms, having two atoms in common.')
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
    identifiers: Optional[Identifiers] = Field(
        None, 
        description = "An optional dictionary of additional identifiers by which this Molecule can be referenced, "
        "such as INCHI, SMILES, SMARTs, etc. See the :class:``Identifiers`` model for more details."
    )
    names: Optional[List[str]] = Field(
        None, 
        description = "A list of atomic label names."
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
    def from_file(cls, filename: Union[FileInput, str], top: Union[FileInput, str] = None, dtype: Optional[str] = None, *, orient: bool = False, **kwargs) -> "Molecule":
        """
        Constructs a molecule object from a file.
        Parameters
        ----------
        filename : str
            The coords filename to build
        top: str
            The topology filename
        dtype : Optional[str], optional
            The type of file to interpret.
        orient : bool, optional
            Orientates the molecule to a standard frame or not.
        **kwargs
            Any additional keywords to pass to the constructor
        Returns
        -------
        Molecule
            A constructed molecule class.
        """
        if not isinstance(filename, FileInput):
            filename = FileInput(path=filename)

        if top and not isinstance(top, FileInput):
            top = FileInput(path=top)
 
        if not dtype:
            if filename.ext in TkMoleculeReaderComponent._extension_maps['qcelem']:
                dtype = TkMoleculeReaderComponent._extension_maps['qcelem'][filename.ext]
                return qcelemental.models.molecule.Molecule.from_file(filename.abs_path, dtype, orient=orient, **kwargs)
        
        if top:
            mol_input = MoleculeReaderInput(file=filename, top_file=top)
        else:
            mol_input = MoleculeReaderInput(file=filename)

        mol = TkMoleculeReaderComponent.compute(mol_input)

        return cls.from_data(mol, dtype=mol.dtype)
        
    @classmethod
    def from_data(cls, data: Any, dtype: Optional[str] = None, *,
        orient: bool = False, validate: bool = None, **kwargs: Dict[str, Any]) -> "Molecule":
        """
        Constructs a molecule object from a data structure.
        Parameters
        ----------
        data: Any
            Data to construct Molecule from
        dtype: Optional[str], optional
            How to interpret the data, if not passed attempts to discover this based on input type.
        orient: bool, optional
            Orientates the molecule to a standard frame or not.
        validate: bool, optional
            Validates the molecule or not.
        **kwargs: Dict[str, Any]
            Additional kwargs to pass to the constructors. kwargs take precedence over data.
        Returns
        -------
        Molecule
            A constructed molecule class.
        """
        if isinstance(data, str):
            try:
                code = ChemCode(code=data)
                mol_input = MoleculeReaderInput(code=code, args={'validate': validate, 'orient': orient, 'kwargs': kwargs})
            except:
                raise ValueError
        elif isinstance(data, ChemCode):
            mol_input = MoleculeReaderInput(code=data, args={'validate': validate, 'orient': orient, 'kwargs': kwargs})
        else:
            # Let's hope this is a toolkit-specific molecule and pass it as data
            mol_input = MoleculeReaderInput(data=data, args={'validate': validate, 'orient': orient, 'kwargs': kwargs})
        
        return MoleculeReaderComponent.compute(mol_input)

    def to_file(self, filename: str, dtype: Optional[str] = None) -> None:
        """Writes the Molecule to a file.
        Parameters
        ----------
        filename : str
            The filename to write to
        dtype : Optional[str], optional
            The type of file to write, attempts to infer dtype from the filename if not provided.
        """
        if not dtype:
            ext = Path(filename).suffix
            for map_name in TkMoleculeReaderComponent._extension_maps:
                if ext in TkMoleculeReaderComponent._extension_maps[map_name]:
                    toolkit = map_name
                    dtype = TkMoleculeReaderComponent._extension_maps[map_name][ext]
                    break
        else:
            for map_name in TkMoleculeReaderComponent._extension_maps:
                if dtype in TkMoleculeReaderComponent._extension_maps[map_name]:
                    toolkit = map_name
                    break

        if toolkit == 'qcelem': 
            super().to_file(filename, dtype)
        elif toolkit == 'rdkit':
            from mmelemental.components.rdkit_component import MoleculeToRDKit
            from rdkit import Chem
            
            rdkmol = MoleculeToRDKit.compute(self)

            if dtype == 'pdb':
                writer = Chem.PDBWriter(filename)
            elif dtype == 'sdf':
                writer = Chem.SDWriter(filename)
            elif dtype == 'smiles':
                writer = Chem.SmilesWriter(filename)
            else:
                raise NotImplementedError(f'File format {dtype} not supported by rdkit.')

            writer.write(rdkmol.mol)
            writer.close()

        elif toolkit == 'parmed':
            from mmelemental.components.parmed_component import MoleculeToParmed
            pmol = MoleculeToParmed.compute(self)
            pmol.mol.save(filename)
        else:
            raise ValueError(f'Data type not yet supported: {dtype}')

    def to_data(self, dtype: str):
        """ Converts Molecule to toolkit-specific molecule (e.g. rdkit). """

        if dtype == 'rdkit':
            from mmelemental.components.rdkit_component import MoleculeToRDKit
            return MoleculeToRDKit.compute(self).mol
        else:
            raise NotImplementedError(f'Data type {dtype} not available.')
