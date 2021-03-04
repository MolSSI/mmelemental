"""
Simulation test for the mmelemental package.
"""
# import pytest
import mmelemental
from mmelemental.models.molecule.mm_mol import Molecule
from mmelemental.models.app.dynamics import DynamicsInput
from mmelemental.models.app.sim_writer import SimWriterInput
from mmelemental.components.io.simwriter_component import SimWriterComponent
import os
from .data import data_dir


def test_mmelemental_md():
    mol = Molecule.from_file(filename=os.path.join(data_dir, "alanine.json"))
    md = DynamicsInput(mol={"mol": mol}, temp=300, press=1.1, temp_method="berendsen")

    sim_input = SimWriterInput(
        model=md, engine=("NAMD", "X.X.X"), filename="input.namd"
    )
    # file = SimWriterComponent.compute(sim_input)
    # Need PSF writer i.e. parmed maybe to get this working
