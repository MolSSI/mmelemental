"""
Simulation test for the mmelemental package.
"""
# import pytest
import mmelemental
from mmelemental.models.molecule.mm_mol import Molecule
from mmelemental.models.proc import ProcInput
import os
from .data import data_mol_dir


def test_mmelemental_md():
    protein = Molecule.from_file(filename=os.path.join(data_mol_dir, "alanine.json"))
    solvent = Molecule.from_file(filename=os.path.join(data_mol_dir, "water.json"))
    proc = ProcInput(
        engine="some_engine",
        engine_version="1.0.0",
        component="mmic_pkg",
        schema_version=0,
    )

    # file = SimWriterComponent.compute(sim_input)
    # Need PSF writer i.e. parmed maybe to get this working
