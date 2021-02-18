"""
Unit and regression test for the mmelemental package.
"""
import pytest
import sys
import os
import parmed
import mmelemental
from mmelemental.models.util.input import FileInput
from mmelemental.models.molecule.mm_mol import Molecule
from mmelemental.models.chem.codes import ChemCode
from mmelemental.models.molecule.io_mol import MolInput, MolOutput
from mmelemental.components.io.constructor_component import (
    MolConstructorComponent,
    ForceFieldConstructorComponent,
)


def test_mmelemental_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmelemental" in sys.modules


def test_mmelemental_moldata():
    groFile = "mmelemental/data/molecules/alanine.gro"
    topFile = "mmelemental/data/molecules/alanine.top"

    mm_mol = Molecule.from_file(groFile, top=topFile)
    assert isinstance(mm_mol, Molecule)

def test_mmelemental_moltop():
    groFile = "mmelemental/data/molecules/alanine.gro"
    topFile = "mmelemental/data/molecules/alanine.top"
    # top = parmed.gromacs.GromacsTopologyFile(topFile.path)
    return Molecule.from_file(groFile, top=topFile)


@pytest.mark.skip(reason="Need rdkit installed to handle codes for now.")
def test_mmelemental_codes():
    smiles = ChemCode(code="CCCC")
    inputs = MolInput(code=smiles)
    return MolConstructorComponent.compute(inputs)


def test_mmelemental_molfiles():
    for ext in ["pdb", "gro"]:
        pdbFile = f"mmelemental/data/molecules/alanine.{ext}"

        mol = Molecule.from_file(pdbFile)

        mol.to_file("mol.pdb")
        # mol.to_file("mol.gro") -> broken in mmic_parmed, why?!
        # mol.to_file("rdkit.xyz")
        # mol.to_file('rdkit.smiles')

        os.remove("mol.pdb")
        # os.remove("mol.gro")
        # os.remove("rdkit.xyz")
        # os.remove('rdkit.smiles')

    return mol
