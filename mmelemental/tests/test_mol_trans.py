"""
Unit and regression test for the mmelemental package.
"""
import pytest
import sys
import os
import mmelemental
from mmelemental.models.util.input import FileInput
from mmelemental.models.molecule.mm_mol import Molecule
import mm_data

try:
    import mmic_translator

    translators = mmic_translator.components.TransComponent.installed_comps()
except Exception:
    translators = []


def pytest_generate_tests(metafunc):
    if "translator" in metafunc.fixturenames:
        metafunc.parametrize("translator", translators)


def test_mmelemental_moldata(translator):
    groFile = mm_data.mols["alanine.gro"]

    mm_mol = Molecule.from_file(groFile, translator=translator)
    assert isinstance(mm_mol, Molecule)


def test_mmelemental_moltop(translator):
    topFile = mm_data.ffs["alanine.top"]
    groFile = mm_data.mols["alanine.gro"]
    mm_mol = Molecule.from_file(groFile, topFile, translator=translator)
    assert mm_mol.connectivity is not None


def test_mmelemental_mol_tofile(translator):
    for ext in ["pdb", "gro"]:
        pdbFile = mm_data.mols[f"alanine.{ext}"]
        mol = Molecule.from_file(pdbFile, translator=translator)

        mol.to_file("mol.pdb")
        mol.to_file("mol.json")  # , indent=2)
        mol.to_file("mol.gro")

        os.remove("mol.pdb")
        os.remove("mol.json")
        os.remove("mol.gro")
    return mol
