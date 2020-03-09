try:
    from rdkit import rdBase, Chem
    from rdkit.Chem import AllChem
except:
    raise ModuleNotFoundError('Make sure rdkit is installed for code validation.')

from qcelemental.models import ProtoModel
from typing import List, Optional, Any, Dict, Tuple

from mmcomponents.components.blueprints.generic_component import GenericComponent

from mmelemental.models.util.input import FileInput
from mmelemental.models.molecule.rdkit_molecule import RDKitMolecule
from mmelemental.models.chem.codes import ChemCode

class MMoleculeReaderInput(ProtoModel):
    file: FileInput = None
    code: ChemCode = None

class MMoleculeReader(GenericComponent):

    _extension_map = {
        ".npy": "numpy",
        ".json": "json",
        ".xyz": "xyz",
        ".psimol": "psi4",
        ".psi4": "psi4",
        ".msgpack": "msgpack",
        ".pdb": "pdb"
    }

    @classmethod
    def input(cls):
        return MMoleculeReaderInput

    @classmethod
    def output(cls):
        return RDKitMolecule

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None) -> Tuple[bool, Dict[str, Any]]:
        
        if isinstance(inputs, dict):
            inputs = MMoleculeReader.input()(**inputs)

        rdkitMol = self.parse_input(inputs)

        return True, RDKitMolecule(mol=rdkitMol) 

    def parse_input(self, inputs: Dict[str, Any], **args) -> Chem.rdchem.Mol:
        """ Creates an instance of rdkit.Chem.Mol by parsing an input file (pdb, etc.) or a chemical code (smiles, etc.)
        """

        # Construct RDKit molecule from a file
        if inputs.file:
            ext = MMoleculeReader._extension_map[inputs.file.ext]
            filename = inputs.file.path

            if ext == 'pdb':
                mol = Chem.MolFromPDBFile(filename, sanitize=False, removeHs=False)
            elif ext == 'rdkMol':
                mol = Chem.MolFromMolFile(filename)
            elif ext == 'mol2':
                mol = Chem.MolFromMol2File(filename)
            elif ext == 'tpl':
                mol = Chem.MolFromTPLFile(filename)
            elif ext == 'sdf':
                mols = Chem.SDMolSupplier(filename)

                if len(mols) > 1:
                    raise ValueError("SDF file should contain a single molecule")
                else:
                    mol = mols[0] # should we support multiple molecules?

            else:
                raise ValueError(f"Unrecognized file type: {ext}")

        # construct RDKit molecule from identifiers
        elif inputs.code:
            codeType = inputs.code.codeType
            function = getattr(Chem, f"MolFrom{codeType}") # should work since validation already done by ChemCode!
            mol = function(inputs.code.code)
            mol = self._gen3D(mol)
        else:
            raise ValueError('Missing input file or code.')

        # Do cleanup if requested
        # TODO: must check for existence of residues
        if args.get('removeResidues'):
            mol = self._removeResidues(args['removeResidues'])

        return mol

    def _gen3D(self, mol, nConformers=1):
        """ Generates 3D coords for a 2D molecule. Should be called only when instantiating a Molecule object.

        :note: a single unique molecule is assumed. 
        """
        mol = Chem.AddHs(mol)
        # create n conformers for molecule
        confargs = AllChem.EmbedMultipleConfs(mol, nConformers)

        # Energy optimize
        for confId in confargs:
            AllChem.UFFOptimizeMolecule(mol, confId=confId)

        return mol

    def _removeResidues(self, residues: List[str]) -> Chem.rdchem.Mol:
        atoms = self._mol.GetAtoms()
        RWmol = Chem.RWMol(mol)

        for atom in atoms:
            if atom.GetPDBResidueInfo().GetResidueName() in residues:
                RWmol.RemoveAtom(atom.GetIdx())

        return Chem.Mol(RWmol)
