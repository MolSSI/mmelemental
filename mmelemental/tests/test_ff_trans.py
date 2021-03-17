"""
Forcefield tests for the mmelemental package.
"""
import pytest
from mmelemental.models.forcefield import ForceField
import sys
import os
from .data import data_ff_dir

try:
    import mmic_translator

    translators = mmic_translator.components.TransComponent.installed()
except Exception:
    translators = []


def pytest_generate_tests(metafunc):
    if "translator" in metafunc.fixturenames:
        metafunc.parametrize("translator", translators)


def test_mmelemental_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmelemental" in sys.modules


def test_mmelemental_moldata(translator):
    topFile = os.path.join(data_ff_dir, "alanine.top")

    mm_ff = ForceField.from_file(topFile, translator=translator)
    assert isinstance(mm_ff, ForceField)
