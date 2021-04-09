import pytest
from mmelemental.models.chem.codes import ChemCode
from mmelemental.models.molecule.mm_mol import Molecule
import os
from mm_data import data_mol_dir


@pytest.mark.skip(reason="Need rdkit installed to handle codes for now.")
def test_mmelemental_codes():
    smiles = ChemCode(code="CCCC")
    inputs = MolInput(code=smiles)
    return MolConstructorComponent.compute(inputs)


def test_mmelemental_json():
    jsonFile = os.path.join(data_mol_dir, "alanine.json")
    mm_mol = Molecule.from_file(jsonFile)
    assert isinstance(mm_mol, Molecule)
