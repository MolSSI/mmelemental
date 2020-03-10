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

    def gen3D(self, nConformers=1):
        """ Generates 3D coords for a molecule. Should be called only when instantiating a Molecule object.

        :note: a single unique molecule is assumed. 
        """
        rdkmol = Chem.AddHs(self.mol)
        # create n conformers for molecule
        confargs = AllChem.EmbedMultipleConfs(rdkmol, nConformers)

        # Energy optimize
        for confId in confargs:
            AllChem.UFFOptimizeMolecule(rdkmol, confId=confId)

        return rdkmol

    def remove_residues(self, residues: List[str]) -> Chem.rdchem.Mol:
        atoms = self._mol.GetAtoms()
        RWmol = Chem.RWMol(self.mol)

        for atom in atoms:
            if atom.GetPDBResidueInfo().GetResidueName() in residues:
                RWmol.RemoveAtom(atom.GetIdx())

        return Chem.Mol(RWmol)

    @classmethod
    def build_mol(cls, inputs: Dict[str, Any], dtype: str) -> Chem.rdchem.Mol:
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
            codeType = inputs.code.codeType
            function = getattr(Chem, f"MolFrom{codeType}") # should work since validation already done by ChemCode!
            rdkmol = function(inputs.code.code)
        else:
            raise ValueError('Missing input file or code.')

        rdkmol = RDKitMolecule(mol=rdkmol)

        if inputs.code:
            rdkmol = rdkmol.gen3D()

        # Do cleanup if requested
        # TODO: must check for existence of residues
        if inputs.resremove:
            rdkmol = rdkmol.removeResidues(inputs.resremove)

        return rdkmol

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