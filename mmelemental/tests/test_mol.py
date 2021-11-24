import pytest
from mmelemental.models import Molecule
from cmselemental.util import yaml_import, which_import
import mm_data
import os

using_yaml = pytest.mark.skipif(
    yaml_import() is False,
    reason="Not detecting module pyyaml or ruamel.yaml. Install package if necessary and add to envvar PYTHONPATH",
)

using_nglview = pytest.mark.skipif(
    which_import("nglview", return_bool=True) is False,
    reason="Not detecting module nglview. Install package if necessary and add to envvar PYTHONPATH",
)

serialize_extensions = [
    "json",
    pytest.param("yaml", marks=using_yaml),
]


@pytest.mark.skip(reason="Need rdkit installed to handle codes for now.")
def test_mmelemental_codes():
    smiles = ChemCode(code="CCCC")
    inputs = MolInput(code=smiles)
    return MolConstructorComponent.compute(inputs)


@pytest.mark.parametrize("encoding", serialize_extensions)
def test_mmelemental_serial(encoding):
    file = mm_data.mols[f"alanine.{encoding}"]
    mm_mol = Molecule.from_file(file)
    assert isinstance(mm_mol, Molecule)

    mm_mol.to_file(f"tmp.{encoding}")
    assert os.path.isfile(f"tmp.{encoding}")
    os.remove(f"tmp.{encoding}")


@pytest.mark.parametrize("encoding", serialize_extensions)
def test_mmelemental_hash(encoding):
    file = mm_data.mols[f"alanine.{encoding}"]
    mm_mol = Molecule.from_file(file)
    data = mm_mol.dict()
    mm_mol = Molecule(**data)
    assert isinstance(mm_mol, Molecule)
