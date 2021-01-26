"""
Forcefield tests for the mmelemental package.
"""
import pytest
from mmelemental.models.forcefield.io_ff import FFInput
import sys


def test_mmelemental_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmelemental" in sys.modules

    ff_in = FFInput(file="mmelemental/data/molecules/alanine.top")
    # ff = ParmedToFFComponent.compute(ff_in)
