"""
Unit and regression test for the mmelemental package.
"""
import pytest
import sys
import os
import parmed
import mmelemental
from mmelemental.models.util.input import FileInput
from mmelemental.models.molecule.mm_molecule import Molecule, MolReaderComponent, MolWriterComponent
from mmelemental.models.chem.codes import ChemCode
from mmelemental.models.molecule.io_molecule import MolInput, MolOutput
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

def test_mmelemental_molfiles(debug=True):
    for ext in ['pdb','gro']:
        pdbFile = FileInput(path=f'mmelemental/data/molecules/dialanine.{ext}')

        mol = Molecule.from_file(filename=pdbFile.path)

        if False:
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
        #mol.to_file('rdkit.smiles')
        
        if not debug:
            os.remove('rdkit.pdb')
            #os.remove('rdkit.gro')
            os.remove('rdkit.xyz')
            #os.remove('rdkit.smiles')

    return mol

def test_mmelemental_molio():
    inp = MolInput(file='mmelemental/data/molecules/dialanine.pdb')
    mol = MolReaderComponent.compute(inp)
    out = MolOutput(file='dialanine.pdb', mol=mol)
    fo = MolWriterComponent.compute(out)

    return fo