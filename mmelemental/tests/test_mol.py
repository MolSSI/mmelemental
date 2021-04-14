import pytest
from mmelemental.models.chem.codes import ChemCode
from mmelemental.models.molecule.mm_mol import Molecule
import mm_data


@pytest.mark.skip(reason="Need rdkit installed to handle codes for now.")
def test_mmelemental_codes():
    smiles = ChemCode(code="CCCC")
    inputs = MolInput(code=smiles)
    return MolConstructorComponent.compute(inputs)


def test_mmelemental_json():
    jsonFile = mm_data.mols["alanine.json"]
    mm_mol = Molecule.from_file(jsonFile)
    assert isinstance(mm_mol, Molecule)
