import qcelemental
from qcelemental.models.types import Array
from typing import List, Tuple, Optional, Any, Dict, Union
import os, sys
import random
import string
import numpy
from pydantic import validator, Field, ValidationError
from mmelemental.components.molreader_component import MMoleculeReader
from mmelemental.models.molecule.molreader import MMoleculeReaderInput
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

class MMolecule(qcelemental.models.Molecule):
    """
    An MMSchema representation of a Molecule based on QCSchema. This model contains data for symbols, geometry, 
    connectivity, charges, residues, etc. while also supporting a wide array of I/O and manipulation capabilities.
    Molecule objects geometry, masses, and charges are truncated to 8, 6, and 4 decimal places respectively 
    to assist with duplicate detection.
    """
    symbols: Array[str] = Field(
        None,
        description = "An ordered (nat,) array-like object of atomic elemental symbols of shape (nat,). The index of "
        "this attribute sets atomic order for all other per-atom setting like ``real`` and the first "
        "dimension of ``geometry``. Ghost/Virtual atoms must have an entry in this array-like and are "
        "indicated by the matching the 0-indexed indices in ``real`` field.",
    )
    geometry: Array[float] = Field( 
        None,
        description = "An ordered (nat,3) array-like for XYZ atomic coordinates [a0]. "
        "Atom ordering is fixed; that is, a consumer who shuffles atoms must not reattach the input "
        "(pre-shuffling) molecule schema instance to any output (post-shuffling) per-atom results "
        "(e.g., gradient). Index of the first dimension matches the 0-indexed indices of all other "
        "per-atom settings like ``symbols`` and ``real``."
        "\n"
        "Can also accept array-likes which can be mapped to (nat,3) such as a 1-D list of length 3*nat, "
        "or the serialized version of the array in (3*nat,) shape; all forms will be reshaped to "
        "(nat,3) for this attribute.",
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
        "of a polymer. Order follows atomic indices from 0 till Natoms-1. "
        "\n"
        "E.g. ('ALA', 1) means atom 0 belongs to aminoacid alanine with residue number 1. Residue number >= 1."
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
        description = "An optional dictionary of additional identifiers by which this MMolecule can be referenced, "
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
    def from_file(cls, filename: str, dtype: Optional[str] = None, *, orient: bool = False, **kwargs) -> "MMolecule":
        """
        Constructs a molecule object from a file.
        Parameters
        ----------
        filename : str
            The filename to build
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
        ext = Path(filename).suffix

        if not dtype:
            if ext in MMoleculeReader._extension_maps['qcelem']:
                dtype = MMoleculeReader._extension_maps['qcelem'][ext]
                return qcelemental.models.molecule.Molecule.from_file(filename, dtype, orient=orient, **kwargs)
        
        molinput = MMoleculeReaderInput(file=FileInput(path=filename))
        mol = MMoleculeReader.compute(molinput)

        return cls.from_data(mol, dtype=mol.obj_type)
        
    @classmethod
    def from_data(cls, data: Union[str, Dict[str, Any], numpy.array, bytes], dtype: Optional[str] = None, *,
        orient: bool = False, validate: bool = None, **kwargs: Dict[str, Any]) -> "MMolecule":
        """
        Constructs a molecule object from a data structure.
        Parameters
        ----------
        data: Union[str, Dict[str, Any], numpy.array]
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

        if not dtype:
            try:
                dtype = data.objType
            except:
                raise ValueError('Input data type (dtype) must be specified for class method: from_data.')

        if dtype == "rdkit":
            try:
                from rdkit import Chem
            except:
                raise ModuleNotFoundError('No installation of rdkit found. '
                        'Make sure rdkit is properly installed on your system.')

            from mmelemental.models.molecule.rdkit_molecule import RDKitMolecule, Bond
            assert isinstance(data, RDKitMolecule)

            symbs = [atom.GetSymbol() for atom in data.mol.GetAtoms()]
            residues = [(atom.GetPDBResidueInfo().GetResidueName(), 
                        atom.GetPDBResidueInfo().GetResidueNumber()) 
                        for atom in data.mol.GetAtoms()]
            names = [atom.GetPDBResidueInfo().GetName() for atom in data.mol.GetAtoms()]

            connectivity = []

            for bond in data.mol.GetBonds():
                bondOrder = Bond.orders.index(bond.GetBondType())
                connectivity.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bondOrder))

            geo = data.mol.GetConformer(0).GetPositions()

            input_dict = {'symbols': symbs, 
                          'geometry': geo, 
                          'residues': residues, 
                          'connectivity': connectivity,
                          'names': names}
        else:
            return qcelemental.models.molecule.Molecule.from_data(data, dtype, orient=orient, validate=validate, **kwargs)

        input_dict.update(kwargs)

        return cls(orient=orient, validate=validate, **input_dict)

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
            for map_name in MMoleculeReader._extension_maps:
                if ext in MMoleculeReader._extension_maps[map_name]:
                    toolkit = map_name
                    dtype = MMoleculeReader._extension_maps[map_name][ext]
                    break
        else:
            for map_name in MMoleculeReader._extension_maps:
                if dtype in MMoleculeReader._extension_maps[map_name]:
                    toolkit = map_name
                    break

        if toolkit == 'qcelem': 
            super().to_file(filename, dtype)
        elif toolkit == 'rdkit':
            try:
                from rdkit import Chem
                from mmelemental.models.molecule.rdkit_molecule import MMToRDKit
            except:
                raise ModuleNotFoundError('Make sure rdkit is installed for code validation.')
            
            rdkmol = MMToRDKit.convert(self)

            if dtype == 'pdb':
                writer = Chem.PDBWriter(filename)
            elif dtype == 'something':
                pass
            else:
                raise NotImplementedError(f'File format {dtype} not supported by rdkit.')
            writer.write(rdkmol)
            writer.close()
        else:
            raise ValueError('Data type {dtype} not supported.')