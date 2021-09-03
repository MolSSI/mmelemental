from pydantic import Field
from typing import List, Dict, Any
from .gen_mol import ToolkitMol
from mmelemental.util.decorators import require

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:  # pragma: no cover
    Chem = AllChem = None  # pragma: no cover


class Bond:
    """RDKit-based bond order: {0: unspecified, 1: single, etc., up to 21}"""

    orders = list(Chem.BondType.values.values())


class RDKitMol(ToolkitMol):
    mol: Chem.rdchem.Mol = Field(..., description="Rdkit molecule object.")

    @require("rdkit")
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @property
    def dtype(self):
        return "rdkit"

    @classmethod
    def gen3D(cls, mol, nConformers=1) -> Chem.rdchem.Mol:
        """Generates 3D coords for a molecule. Should be called only when instantiating a Molecule object.

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
    def build(cls, inputs: Dict[str, Any], dtype: str = None) -> "RDKitMol":
        """Creates an instance of RDKitMol object storing rdkit.Chem.Mol.
        This is done by parsing an input file (pdb, ...) or a chemical code (smiles, ...).
        """
        if inputs.file:
            if not dtype:
                dtype = filename.ext
            filename = inputs.file.path
            if dtype == ".pdb" or dtype == "pdb":
                rdkmol = Chem.MolFromPDBFile(filename, sanitize=False, removeHs=False)
            elif dtype == ".mol" or dtype == "mol":
                rdkmol = Chem.MolFromMolFile(filename, sanitize=False, removeHs=False)
            elif dtype == ".mol2" or dtype == "mol2":
                rdkmol = Chem.MolFromMol2File(filename, sanitize=False, removeHs=False)
            elif dtype == ".tpl" or dtype == "tpl":
                rdkmol = Chem.MolFromTPLFile(filename, sanitize=False, removeHs=False)
            elif dtype == ".sdf" or dtype == "sdf":
                rdkmols = Chem.SDMolSupplier(filename, sanitize=False, removeHs=False)

                if len(rdkmols) > 1:
                    raise ValueError("SDF file should contain a single molecule")
                else:
                    rdkmol = rdkmols[0]  # should we support multiple molecules?
            else:
                raise ValueError(f"Unrecognized file type: {dtype}")

        # construct RDKit molecule from identifiers
        elif inputs.code:
            code_type = inputs.code.code_type
            function = getattr(
                Chem, f"MolFrom{code_type}"
            )  # should work since validation already done by ChemCode!
            rdkmol = function(inputs.code.code)
        else:
            raise ValueError("Missing input file or code.")

        if inputs.code:
            rdkmol = cls.gen3D(rdkmol)

        return cls(mol=rdkmol)
