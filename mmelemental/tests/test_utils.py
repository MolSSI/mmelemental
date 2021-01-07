import sys
import os
debug = True

# import models
from mmelemental.models.util.input import FileInput, OpenBabelInput, GrepInput
from mmelemental.models.chem.codes import ChemCode

# Test input file
receptor = FileInput(path=os.path.abspath('mmelemental/data/molecules/dialanine.pdb'))

def test_obabel():
   from mmelemental.components.util.openbabel_component import OpenBabel

   obabel_input = OpenBabelInput(fileInput=receptor, outputExt='pdbqt')
   obabel_output = OpenBabel.compute(obabel_input)

def test_grep():
   from mmelemental.components.util.grep_component import Grep

   grep_input = GrepInput(fileInput=receptor, pattern='ATOM')
   grep_output = Grep.compute(grep_input)
