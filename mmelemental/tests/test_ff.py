"""
Forcefield tests for the mmelemental package.
"""
import pytest
from mmelemental.models.forcefield import ForceField
import sys
import os

data_dir = os.path.join("mmelemental", "data", "molecules")


def test_mmelemental_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmelemental" in sys.modules


def test_mmelemental_read_gmx():
    groFile = os.path.join(data_dir, "alanine.top")

    mm_ff = ForceField.from_file(groFile)
    assert isinstance(mm_ff, ForceField)
