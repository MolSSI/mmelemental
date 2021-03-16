"""
Simulation test for the mmelemental package.
"""
# import pytest
import mmelemental
from mmelemental.models.molecule.mm_mol import Molecule
from mmelemental.models.proc.dynamics import DynamicsInput
import os
from .data import data_dir


def test_mmelemental_md():
    protein = Molecule.from_file(filename=os.path.join(data_dir, "alanine.json"))
    solvent = Molecule.from_file(filename=os.path.join(data_dir, "water.json"))
    md = DynamicsInput(
        mol={"protein": protein, "solvent": solvent},
        forcefield={"protein": "charmm27", "solvent": "spc"},
        temp=300,
        press=1.1,
        temp_method="berendsen",
    )

    # file = SimWriterComponent.compute(sim_input)
    # Need PSF writer i.e. parmed maybe to get this working
