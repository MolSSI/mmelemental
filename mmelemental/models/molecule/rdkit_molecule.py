from qcelemental import models
from pydantic import Field
from typing import List, Dict, Any
from .gen_molecule import ToolkitMolecule

try:
    from rdkit import rdBase, Chem
    from rdkit.Chem import AllChem
except:
    raise ModuleNotFoundError('Make sure rdkit is installed for code validation.')

class Bond:
    """ RDKit-based bond order: {0: unspecified, 1: single, etc., up to 21} """  
    orders = list(Chem.BondType.values.values())

class RDKitMolecule(ToolkitMolecule):
    mol: Chem.rdchem.Mol = Field(None, description = 'Rdkit molecule object.')
        
    @property
    def obj_type(self):
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
    def build_mol(cls, inputs: Dict[str, Any], dtype: str) -> "RDKitMolecule":
        """ Creates an instance of rdkit.Chem.Mol by parsing an input file (pdb, etc.) or a chemical code (smiles, etc.)
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

        # Do cleanup if requested
        # TODO: must check for existence of residues
        if inputs.resremove:
            rdkmol = cls.removeResidues(rdkmol, inputs.resremove)

        return RDKitMolecule(mol=rdkmol)

class MMToRDKit:
    @staticmethod
    def convert(mmol: "MMolecule") -> RDKitMolecule:
        rdkmol = Chem.Mol()
        erdkmol = Chem.EditableMol(rdkmol)

        for index, symb in enumerate(mmol.symbols):
            atom = Chem.Atom(symb)
            resname, resnum = mmol.residues[index]
            name = mmol.names[index]
            residue = Chem.AtomPDBResidueInfo()
            residue.SetResidueName(resname)
            residue.SetName(name)
            residue.SetResidueNumber(resnum)
            residue.SetOccupancy(1.0)
            residue.SetTempFactor(0.0)
            atom.SetMonomerInfo(residue)
            erdkmol.AddAtom(atom)

        for i,j,btype in mmol.connectivity:
            erdkmol.AddBond(i, j, Bond.orders[int(btype)])            

        newmmol = erdkmol.GetMol()
        conf = Chem.Conformer(len(mmol.geometry))
        for i, coords in enumerate(mmol.geometry):
            conf.SetAtomPosition(i, coords)
        newmmol.AddConformer(conf)

        return newmmol 
