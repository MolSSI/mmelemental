import pytest
from pathlib import Path
from mmelemental.models import Molecule
from mmelemental.models.struct import Topology
from cmselemental.util import yaml_import, which_import
import mm_data
import os

using_yaml = pytest.mark.skipif(
    yaml_import() is False,
    reason="Not detecting module pyyaml or ruamel.yaml. Install package if necessary and add to envvar PYTHONPATH",
)

serialize_extensions = [
    "json",
    pytest.param("yaml", marks=using_yaml),
]


@pytest.mark.parametrize("encoding", serialize_extensions)
def test_mmelemental_top_serial(encoding):
    file = mm_data.mols[f"alanine.{encoding}"]
    mm_mol = Molecule.from_file(file)
    assert isinstance(mm_mol, Molecule)

    top = mm_mol.get_topology()
    mol_from_top = Molecule(**top.dict(exclude={"schema_name"}))

    path = Path(f"tmp.{encoding}")
    path.write_text(top.json())
    assert path.is_file()
    path.unlink()


def test_mmelemental_top():
    file = mm_data.mols[f"alanine.json"]
    mm_mol = Molecule.from_file(file)
    assert isinstance(mm_mol, Molecule)

    top_from_mol = mm_mol.get_topology()
    top = Topology(**top_from_mol.dict(exclude={"schema_name"}))

    blob = top.json()
