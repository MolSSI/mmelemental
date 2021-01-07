"""
Unit and regression test for the mmelemental package.
"""
import pytest
import sys
import os
import parmed
import mmelemental
from mmelemental.models.util.input import FileInput
from mmelemental.models.molecule.mm_molecule import Molecule
from mmelemental.models.chem.codes import ChemCode
from mmelemental.models.molecule.mol_reader import MolInput

from mmelemental.components.io.molreader_component import MolReaderComponent
from mmelemental.components.io.constructor_component import MolConstructorComponent, ForceFieldConstructorComponent

def test_mmelemental_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmelemental" in sys.modules

def test_mmelemental_moltop():
    groFile = FileInput(path='mmelemental/data/molecules/dialanine.gro')
    topFile = FileInput(path='mmelemental/data/molecules/dialanine.top')
    #top = parmed.gromacs.GromacsTopologyFile(topFile.path)
    mol = Molecule.from_file(filename=groFile, top=topFile)

@pytest.mark.skip(reason="Need rdkit installed to handle codes for now.")
def test_mmelemental_codes():
    smiles = ChemCode(code='CCCC')
    inputs = MolInput(code=smiles)
    mol = MolConstructorComponent.compute(inputs)

def test_mmelemental_molfiles(debug=False):
    for ext in ['pdb']:
        pdbFile = FileInput(path=f'mmelemental/data/molecules/dialanine.{ext}')

        mol = Molecule.from_file(filename=pdbFile.path)

        if debug:
            print("Molecule info:")
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

        mol.to_file('rdkit.pdb')
        #mol.to_file('rdkit.gro')
        mol.to_file('rdkit.xyz')
        mol.to_file('rdkit.smiles')
        
        if not debug:
            os.remove('rdkit.pdb')
            #os.remove('rdkit.gro')
            os.remove('rdkit.xyz')
            os.remove('rdkit.smiles')

    return mol