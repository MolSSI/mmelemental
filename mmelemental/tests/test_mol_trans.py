"""
Unit and regression test for the mmelemental package.
"""
import pytest
import sys
import os
import mmelemental
from mmelemental.models.util.input import FileInput
from mmelemental.models.molecule.mm_mol import Molecule

from mmelemental.models.molecule.io_mol import MolInput, MolOutput
from mmelemental.components.io.constructor_component import (
    MolConstructorComponent,
    ForceFieldConstructorComponent,
)

from .data import data_mol_dir, data_ff_dir

try:
    import mmic_translator

    translators = mmic_translator.components.TransComponent.installed()
except Exception:
    translators = []


def pytest_generate_tests(metafunc):
    if "translator" in metafunc.fixturenames:
        metafunc.parametrize("translator", translators)


def test_mmelemental_moldata(translator):
    groFile = os.path.join(data_mol_dir, "alanine.gro")

    mm_mol = Molecule.from_file(groFile, translator=translator)
    assert isinstance(mm_mol, Molecule)


def test_mmelemental_moltop(translator):
    topFile = os.path.join(data_ff_dir, "alanine.top")
    groFile = os.path.join(data_mol_dir, "alanine.gro")
    mm_mol = Molecule.from_file(groFile, topFile, translator=translator)
    assert mm_mol.connectivity is not None


def test_mmelemental_mol_tofile(translator):
    for ext in ["pdb", "gro"]:
        pdbFile = os.path.join(data_mol_dir, f"alanine.{ext}")

        mol = Molecule.from_file(pdbFile, translator=translator)

        mol.to_file("mol.pdb")
        mol.to_file("mol.json")  # , indent=2)
        # mol.to_file("mol.gro") -> broken in mmic_parmed, why?!
        # mol.to_file("rdkit.xyz")
        # mol.to_file('rdkit.smiles')

        os.remove("mol.pdb")
        os.remove("mol.json")
        # os.remove("mol.gro")
        # os.remove("rdkit.xyz")
        # os.remove('rdkit.smiles')

    return mol
