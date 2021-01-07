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
from mmelemental.models.input.dynamics import DynamicsInput
from mmelemental.components.io.simwriter_component import SimWriterComponent
from mmelemental.models.input.sim_writer import SimWriterInput

def test_mmelemental_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmelemental" in sys.modules

def test_mmelemental_molgro():
    groFile = FileInput(path='mmelemental/data/molecules/dialanine.gro')
    topFile = FileInput(path='mmelemental/data/molecules/dialanine.top')

    # top = parmed.gromacs.GromacsTopologyFile(topFile.path)

    return Molecule.from_file(filename=groFile.path, top_file=topFile)

mol = test_mmelemental_molgro()
md = DynamicsInput(
		mol   = mol,
		temp  = 300,
		press = 1.1,
                temp_method = 'berendsen'
	)

print(md) 

sim_input = SimWriterInput(model=md, engine=('NAMD', 'X.X.X'), filename='input.namd')
file = SimWriterComponent.compute(sim_input)
