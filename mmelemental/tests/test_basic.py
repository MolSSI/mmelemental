"""
Basic sanity test for the mmelemental package.
"""
# import pytest
import mmelemental
from mmelemental.models.util.output import FileOutput
from mmelemental.models.molecule import Molecule


def test_mmelemental_imported():
    """Sample test, will always pass so long as import statement worked"""
    import sys

    assert "mmelemental" in sys.modules


def test_mmelemental_mol_schema():
    import os

    mol = Molecule.from_file("mmelemental/data/molecules/alanine.json")
    with FileOutput(path="mol.json", clean=True) as fo:
        mol.to_file(fo.abs_path)
