"""
Unit and regression test for the mmelemental package.
"""
#import pytest
import sys
import os
import parmed
import mmelemental
from mmelemental.models.util.input import FileInput
from mmelemental.models.molecule.mm_molecule import MMolecule
from mmelemental.models.chem.codes import ChemCode
from mmelemental.models.molecule.mol_reader import MMoleculeReaderInput

from mmelemental.components.molreader_component import MMoleculeReaderComponent
from mmelemental.components.constructor_component import MolConstructorComponent, ForceFieldConstructorComponent

def test_mmelemental_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmelemental" in sys.modules

def test_mmelemental_molgro(debug=True):
    groFile = FileInput(path=os.path.abspath('mmelemental/data/molecules/dialanine.gro'))
    topFile = FileInput(path=os.path.abspath('mmelemental/data/molecules/dialanine.top'))
    
    top = parmed.gromacs.GromacsTopologyFile(topFile.path)

    return MMolecule.from_file(filename=groFile.path)

def test_mmelemental_molpdb(debug=True):
    pdbFile = FileInput(path=os.path.abspath('mmelemental/data/molecules/dialanine.pdb'))
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
    mol.to_file('tmp.xyz')
    mol.to_file('tmp.smiles')
    os.remove('tmp.pdb')
    os.remove('tmp.xyz')
    os.remove('tmp.smiles')

def test_mmelemental_component():

    smiles = ChemCode(code='CCCC')
    mol = MolConstructorComponent.compute(MMoleculeReaderInput(code=smiles))
    mol = MMoleculeReaderComponent.compute(MMoleculeReaderInput(code=smiles))

test_mmelemental_imported()
test_mmelemental_molpdb()
test_mmelemental_component()
