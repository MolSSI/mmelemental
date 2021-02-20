"""
Forcefield tests for the mmelemental package.
"""
import pytest
from mmelemental.models.forcefield import ForceField
import sys


def test_mmelemental_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmelemental" in sys.modules
