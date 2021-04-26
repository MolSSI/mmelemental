"""
Simulation test for the mmelemental package.
"""
# import pytest
import mmelemental
from mmelemental.models.molecule.mm_mol import Molecule
from mmelemental.models.proc import ProcInput
import mm_data


def test_mmelemental_md():
    protein = Molecule.from_file(filename=mm_data.mols["alanine.json"])
    solvent = Molecule.from_file(filename=mm_data.mols["water-mol.json"])
    proc = ProcInput(
        engine="some_engine",
        engine_version="1.0.0",
        component="mmic_pkg",
        schema_version=0,
    )

    # file = SimWriterComponent.compute(sim_input)
    # Need PSF writer i.e. parmed maybe to get this working
