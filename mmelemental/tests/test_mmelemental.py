"""
Unit and regression test for the mmelemental package.
"""

# Import package, test suite, and other packages as needed
import mmelemental
import pytest
import sys

def test_mmelemental_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmelemental" in sys.modules
