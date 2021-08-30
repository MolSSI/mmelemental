import pytest
from pathlib import Path
from mmelemental.models.molecule.mm_mol import Molecule
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
