from qcelemental import models
from pydantic import Field, validator
from typing import List, Dict, Any
from .gen_molecule import ToolkitMolecule

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except:
    raise ModuleNotFoundError('Make sure rdkit is installed.')

class Bond:
    """ RDKit-based bond order: {0: unspecified, 1: single, etc., up to 21} """  
    orders = list(Chem.BondType.values.values())

class RDKitMolecule(ToolkitMolecule):
    mol: Chem.rdchem.Mol = Field(..., description = 'Rdkit molecule object.')

    @property
    def dtype(self):
        return 'rdkit'

    @classmethod
    def gen3D(cls, mol, nConformers=1) -> Chem.rdchem.Mol:
        """ Generates 3D coords for a molecule. Should be called only when instantiating a Molecule object.

        :note: a single unique molecule is assumed. 
        """
        rdkmol = Chem.AddHs(mol)
        # create n conformers for molecule
        confargs = AllChem.EmbedMultipleConfs(rdkmol, nConformers)

        # Energy optimize
        for confId in confargs:
            AllChem.UFFOptimizeMolecule(rdkmol, confId=confId)

        return rdkmol

    @classmethod
    def remove_residues(cls, mol, residues: List[str]) -> Chem.rdchem.Mol:
        atoms = mol.GetAtoms()
        RWmol = Chem.RWMol(mol)

        for atom in atoms:
            if atom.GetPDBResidueInfo().GetResidueName() in residues:
                RWmol.RemoveAtom(atom.GetIdx())

        return Chem.Mol(RWmol)

    @classmethod
    def build(cls, inputs: Dict[str, Any], dtype: str) -> "RDKitMolecule":
        """ Creates an instance of RDKitMolecule object storing rdkit.Chem.Mol. 
        This is done by parsing an input file (pdb, ...) or a chemical code (smiles, ...).
        """
        if inputs.file:
            filename = inputs.file.path
            if dtype == 'pdb':
                rdkmol = Chem.MolFromPDBFile(filename, sanitize=False, removeHs=False)
            elif dtype == 'mol':
                rdkmol = Chem.MolFromMolFile(filename, sanitize=False, removeHs=False)
            elif dtype == 'mol2':
                rdkmol = Chem.MolFromMol2File(filename, sanitize=False, removeHs=False)
            elif dtype == 'tpl':
                rdkmol = Chem.MolFromTPLFile(filename, sanitize=False, removeHs=False)
            elif dtype == 'sdf':
                rdkmols = Chem.SDMolSupplier(filename, sanitize=False, removeHs=False)

                if len(rdkmols) > 1:
                    raise ValueError("SDF file should contain a single molecule")
                else:
                    rdkmol = rdkmols[0] # should we support multiple molecules?

            else:
                raise ValueError(f"Unrecognized file type: {ext}")

        # construct RDKit molecule from identifiers
        elif inputs.code:
            code_type = inputs.code.code_type
            function = getattr(Chem, f"MolFrom{code_type}") # should work since validation already done by ChemCode!
            rdkmol = function(inputs.code.code)
        else:
            raise ValueError('Missing input file or code.')

        if inputs.code:
            rdkmol = cls.gen3D(rdkmol)

        return cls(mol=rdkmol)