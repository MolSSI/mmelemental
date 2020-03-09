"""
Unit and regression test for the mmelemental package.
"""
import pytest
import sys
import os

import mmelemental
from mmelemental.models.util.input import FileInput
from mmelemental.models.molecule.mm_molecule import MMolecule
from mmelemental.models.chem.codes import ChemCode
from mmelemental.components.mm_reader import MMoleculeReader, MMoleculeReaderInput

def test_mmelemental_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmelemental" in sys.modules

def test_mmelemental_models(debug):
    pdbFile = FileInput(path=os.path.abspath('mmelemental/data/molecules/dialanine.pdb'))
    sdfFile = FileInput(path=os.path.abspath('mmelemental/data/molecules/ligand.sdf'))
    smiles = ChemCode(code='CCC')

    mol = MMolecule.from_file(filename=pdbFile.path)

    if debug:
        print("MMolecule info:")
        print("===============")
        print('\n')

        print("Bonds:")
        print("======")
        print(mol.connectivity)
        print('\n')

        print("Residues:")
        print("==========")
        print(mol.residues)
        print('\n')

        print("Positions:")
        print("==========")
        print(mol.geometry)
        print('\n')

        print("Atom Names:")
        print("==========")
        print(mol.names)
        print('\n')

    mol.to_file('tmp.pdb')

test_mmelemental_imported()
test_mmelemental_models(True)
