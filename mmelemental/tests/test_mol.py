"""
Unit and regression test for the mmelemental package.
"""
import pytest
import sys
import os
import parmed
import mmelemental
from mmelemental.models.util.input import FileInput
from mmelemental.models.molecule.mm_mol import Molecule
from mmelemental.models.chem.codes import ChemCode
from mmelemental.models.molecule.io_mol import MolInput, MolOutput
from mmelemental.components.io.constructor_component import (
    MolConstructorComponent,
    ForceFieldConstructorComponent,
)


translators = mmelemental.components.trans.TransComponent.installed()
data_dir = os.path.join("mmelemental", "data", "molecules")


def pytest_generate_tests(metafunc):
    if "translator" in metafunc.fixturenames:
        metafunc.parametrize("translator", translators)


def test_mmelemental_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmelemental" in sys.modules


def test_mmelemental_moldata(translator):
    groFile = os.path.join(data_dir, "alanine.gro")

    mm_mol = Molecule.from_file(groFile, translator=translator)
    assert isinstance(mm_mol, Molecule)


def test_mmelemental_moltop(translator):
    topFile = os.path.join(data_dir, "alanine.top")
    groFile = os.path.join(data_dir, "alanine.gro")

    import importlib
    from pathlib import Path

    mod = importlib.import_module(translator)

    if Path(topFile).suffix in mod.ffread_ext_maps:
        mm_mol = Molecule.from_file(groFile, topFile, translator=translator)
        assert mm_mol.connectivity is not None


@pytest.mark.skip(reason="Need rdkit installed to handle codes for now.")
def test_mmelemental_codes():
    smiles = ChemCode(code="CCCC")
    inputs = MolInput(code=smiles)
    return MolConstructorComponent.compute(inputs)


def test_mmelemental_json():
    jsonFile = os.path.join(data_dir, "alanine.json")
    mm_mol = Molecule.from_file(jsonFile)
    assert isinstance(mm_mol, Molecule)


def test_mmelemental_mol_tofile(translator):
    for ext in ["pdb", "gro"]:
        pdbFile = os.path.join(data_dir, f"alanine.{ext}")

        mol = Molecule.from_file(pdbFile, translator=translator)

        mol.to_file("mol.pdb")
        mol.to_file("mol.json", indent=2)
        # mol.to_file("mol.gro") -> broken in mmic_parmed, why?!
        # mol.to_file("rdkit.xyz")
        # mol.to_file('rdkit.smiles')

        os.remove("mol.pdb")
        os.remove("mol.json")
        # os.remove("mol.gro")
        # os.remove("rdkit.xyz")
        # os.remove('rdkit.smiles')

    return mol
