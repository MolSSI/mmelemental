"""
Unit and regression test for the mmelemental package.
"""
import pytest
import sys
import os
import parmed
import mmelemental
from mmelemental.models.util.input import FileInput
from mmelemental.models.molecule.mm_mol import (
    Mol,
)
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
    groFile = FileInput(path="mmelemental/data/molecules/alanine.gro")
    topFile = FileInput(path="mmelemental/data/molecules/alanine.top")

    mm_mol = Mol.from_file(groFile, top=topFile)
    assert(isinstance(mm_mol, Mol))

    mda_mol = mm_mol.to_data(dtype="MDAnalysis")
    assert(isinstance(mda_mol.mol, mda_mol.dtype))

def test_mmelemental_moltop():
    groFile = FileInput(path="mmelemental/data/molecules/alanine.gro")
    topFile = FileInput(path="mmelemental/data/molecules/alanine.top")
    # top = parmed.gromacs.GromacsTopologyFile(topFile.path)
    mol = Mol.from_file(groFile, top=topFile)


@pytest.mark.skip(reason="Need rdkit installed to handle codes for now.")
def test_mmelemental_codes():
    smiles = ChemCode(code="CCCC")
    inputs = MolInput(code=smiles)
    mol = MolConstructorComponent.compute(inputs)


def test_mmelemental_molfiles(debug=True):
    for ext in ["pdb", "gro"]:
        pdbFile = FileInput(path=f"mmelemental/data/molecules/alanine.{ext}")

        mol = Mol.from_file(pdbFile.path)

        if False:
            print("Molecule info:")
            print("===============")
            print("\n")

            print("Bonds:")
            print("======")
            print(mol.connectivity)
            print("\n")

            print("Residues:")
            print("==========")
            print(mol.residues)
            print("\n")

            print("Positions:")
            print("==========")
            print(mol.geometry)
            print("\n")

            print("Atom Names:")
            print("==========")
            print(mol.names)
            print("\n")

        mol.to_file("rdkit.pdb")
        # mol.to_file('rdkit.gro')
        mol.to_file("rdkit.xyz")
        # mol.to_file('rdkit.smiles')

        if not debug:
            os.remove("rdkit.pdb")
            # os.remove('rdkit.gro')
            os.remove("rdkit.xyz")
            # os.remove('rdkit.smiles')

    return mol