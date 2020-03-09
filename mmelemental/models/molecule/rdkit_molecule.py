from qcelemental import models

try:
    from rdkit import rdBase, Chem
    from rdkit.Chem import AllChem
except:
    raise ModuleNotFoundError('Make sure rdkit is installed for code validation.')

class Bond:
    """ RDKit-based bond order: {0: unspecified, 1: single, etc., up to 21} """  
    orders = list(Chem.BondType.values.values())

class RDKitMolecule(models.ProtoModel):
    mol: Chem.rdchem.Mol = None

    class Config:
        arbitrary_types_allowed = True

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
