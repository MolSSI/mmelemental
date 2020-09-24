from mmcomponents.components.blueprints.generic_component import GenericComponent
from mmelemental.models.util.output import FileOutput
from mmelemental.models.molecule.rdkit_molecule import RDKitMolecule
from mmelemental.models.molecule.mm_molecule import Molecule
from mmelemental.models.molecule.rdkit_molecule import Bond
from typing import Dict, Any, List, Tuple, Optional

try:
    from rdkit import Chem
except:
    raise ModuleNotFoundError('Make sure rdkit is installed.')

class MoleculeToRDKit(GenericComponent):
    """ A model for converting Molecule to RDKIT molecule object. """
    @classmethod
    def input(cls):
        return Molecule

    @classmethod
    def output(cls):
        return RDKitMolecule

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, RDKitMolecule]:

        rdkmol = Chem.Mol()
        erdkmol = Chem.EditableMol(rdkmol)
        mmol = inputs

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

        return True, RDKitMolecule(mol=newmmol)
