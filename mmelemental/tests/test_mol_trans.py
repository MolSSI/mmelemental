"""
Unit and regression test for the mmelemental package.
"""
import pytest
import sys
import os
import importlib
import mmelemental
import mm_data

try:
    import mmic_translator

    translators = mmic_translator.components.TransComponent.installed_comps_model(
        "Molecule"
    )
except Exception:
    translators = []


def pytest_generate_tests(metafunc):
    if "translator" in metafunc.fixturenames:
        metafunc.parametrize("translator", translators)


def test_mmelemental_mol(translator):
    trans_mod = importlib.import_module(translator)
    for ext in ["pdb", "gro", "xyz"]:
        if ext in trans_mod.molread_ext_maps:
            molFile = mm_data.mols[f"alanine.{ext}"]
            mol = mmelemental.models.Molecule.from_file(molFile, translator=translator)

        if ext in trans_mod.molwrite_ext_maps:
            mol.to_file(f"mol.{ext}")
            os.remove(f"mol.{ext}")


def test_mmelemental_moltop(translator):
    trans_mod = importlib.import_module(translator)
    for ext in ["pdb", "gro"]:
        if ext in trans_mod.molread_ext_maps:
            molFile = mm_data.mols[f"alanine.{ext}"]
            topFile = mm_data.ffs["alanine.top"]
            mol = mmelemental.models.Molecule.from_file(
                molFile, topFile, translator=translator
            )
            assert mol.connectivity is not None

        if ext in trans_mod.molwrite_ext_maps:
            mol.to_file(f"mol.{ext}")
            os.remove(f"mol.{ext}")
