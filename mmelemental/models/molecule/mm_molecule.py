import qcelemental
from qcelemental.models.types import Array
from typing import List, Tuple, Optional, Any, Dict, Union
import os, sys
import random
import string
import numpy
from pydantic import validator, Field, ValidationError
from mmelemental.components.io.molreader_component import TkMolReaderComponent, MolReaderComponent
from mmelemental.models.molecule.mol_reader import MolReaderInput
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


import unyt
from typing import Any, Dict
import numpy as np

class TypedArray(unyt.unyt_array):
    @classmethod
    def __get_validators__(cls):
        yield cls.validate

    @classmethod
    def validate(cls, v):
        try:
            v = unyt.unyt_array(v, dtype=cls._dtype)
        except ValueError:
            raise ValueError("Could not cast {} to NumPy Array!".format(v))

        return v

    @classmethod
    def __modify_schema__(cls, field_schema: Dict[str, Any]) -> None:
        dt = cls._dtype
        if dt is int or np.issubdtype(dt, np.integer):
            items = {"type": "number", "multipleOf": 1.0}
        elif dt is float or np.issubdtype(dt, np.floating):
            items = {"type": "number"}
        elif dt is str or np.issubdtype(dt, np.string_):
            items = {"type": "string"}
        elif dt is bool or np.issubdtype(dt, np.bool_):
            items = {"type": "boolean"}
        field_schema.update(type="array", items=items)

class ArrayMeta(type):
    def __getitem__(self, dtype):
        return type("UnytArray", (TypedArray,), {"_dtype": dtype})

class UnytArray(np.ndarray, metaclass=ArrayMeta):
    pass

class Molecule(qcelemental.models.Molecule):
    """
    An MMSchema representation of a Molecule based on QCSchema. This model contains data for symbols, geometry, 
    connectivity, charges, residues, etc. while also supporting a wide array of I/O and manipulation capabilities.
    Molecule objects geometry, masses, and charges are truncated to 8, 6, and 4 decimal places respectively 
    to assist with duplicate detection.
    """
    symbols: Optional[Array[str]] = Field(
        None,
        description = "An ordered (natom,) array-like object of atomic elemental symbols. The index of "
        "this attribute sets atomic order for all other per-atom setting like ``real`` and the first "
        "dimension of ``geometry``. Ghost/Virtual atoms must have an entry in this array-like and are "
        "indicated by the matching the 0-indexed indices in ``real`` field.",
    )
    geometry: Union[Array[float], UnytArray[float]] = Field(
        None,
        description="An ordered (nat,3) array-like for XYZ atomic coordinates [a0]. "
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
            filename = FileInput(path=filename)

        if top and not isinstance(top, FileInput):
            top = FileInput(path=top)
 
        if not dtype:
            if filename.ext in TkMolReaderComponent._extension_maps['qcelement']:
                dtype = TkMolReaderComponent._extension_maps['qcelement'][filename.ext]
                return qcelemental.models.molecule.Molecule.from_file(filename.abs_path, dtype, orient=orient, **kwargs)
        if top:
            mol_input = MolReaderInput(file=filename, top_file=top)
        else:
            mol_input = MolReaderInput(file=filename)

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
                mol_input = MolReaderInput(code=code, args={'validate': validate, 'orient': orient, 'kwargs': kwargs})
            except:
                raise ValueError
        elif isinstance(data, ChemCode):
            mol_input = MolReaderInput(code=data, args={'validate': validate, 'orient': orient, 'kwargs': kwargs})
        else:
            # Let's hope this is a toolkit-specific molecule and pass it as data
            mol_input = MolReaderInput(data=data, args={'validate': validate, 'orient': orient, 'kwargs': kwargs})
        
        return MolReaderComponent.compute(mol_input)

    def to_file(self, filename: str, dtype: Optional[str] = None) -> None:
        """ Writes the Molecule to a file.
        Parameters
        ----------
        filename : str
            The filename to write to
        dtype : Optional[str], optional
            The type of file to write, attempts to infer dtype from the filename if not provided.
        """
        if not dtype:
            ext = Path(filename).suffix
            for map_name in TkMolReaderComponent._extension_maps:
                if ext in TkMolReaderComponent._extension_maps[map_name]:
                    toolkit = map_name
                    dtype = TkMolReaderComponent._extension_maps[map_name][ext]
                    break
        else:
            for map_name in TkMolReaderComponent._extension_maps:
                if dtype in TkMolReaderComponent._extension_maps[map_name]:
                    toolkit = map_name
                    break

        if toolkit == 'qcelement': 
            super().to_file(filename, dtype)
        elif toolkit == 'rdkit':
            from mmelemental.components.trans.rdkit_component import MoleculeToRDKit
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
            from mmelemental.components.trans.parmed_component import MoleculeToParmed
            pmol = MoleculeToParmed.compute(self)
            pmol.mol.save(filename)
        else:
            raise ValueError(f'Data type not yet supported: {dtype}')

    def to_data(self, dtype: str):
        """ Converts Molecule to toolkit-specific molecule (e.g. rdkit). """

        if dtype == 'rdkit':
            from mmelemental.components.trans.rdkit_component import MoleculeToRDKit
            return MoleculeToRDKit.compute(self).mol
        else:
            raise NotImplementedError(f'Data type {dtype} not available.')

    def get_molecular_formula(self, order: str = "alphabetical") -> str:
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

        return super().get_molecular_formula(order)