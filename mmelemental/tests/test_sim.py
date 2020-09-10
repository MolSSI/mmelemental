"""
Simulation test for the mmelemental package.
"""
#import pytest
import sys
import os
import parmed
import mmelemental
from mmelemental.models.util.input import FileInput
from mmelemental.models.molecule.mm_molecule import Molecule
from mmelemental.models.molecule.mol_reader import MoleculeReaderInput
from mmelemental.models.sim.md import MolecularDynamics

def test_mmelemental_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmelemental" in sys.modules

def test_mmelemental_molpdb():
    groFile = FileInput(path=os.path.abspath('mmelemental/data/molecules/dialanine.pdb'))
    topFile = FileInput(path=os.path.abspath('mmelemental/data/molecules/dialanine.top'))

    # top = parmed.gromacs.GromacsTopologyFile(topFile.path)

    return Molecule.from_file(filename=groFile.path)

mol = test_mmelemental_molpdb()
md = MolecularDynamics(mol=mol)

print(md)
