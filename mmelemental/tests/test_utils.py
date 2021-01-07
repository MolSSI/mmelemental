import sys
import os
debug = True

# import models
from mmelemental.models.util.input import FileInput, OpenBabelInput, GrepInput
from mmelemental.models.chem.codes import ChemCode

# Test input file
receptor = FileInput(path=os.path.abspath('mmelemental/data/molecules/dialanine.pdb'))

# Import components
from mmelemental.components.util.openbabel_component import OpenBabel
from mmelemental.components.util.grep_component import Grep

# Test for openbabel
obabel_input = OpenBabelInput(fileInput=receptor, outputExt='pdbqt')
obabel_output = OpenBabel.compute(obabel_input)

if debug:
	print("==============================")
	print("OBABEL OUTPUT:")
	print("==============================")
	print(obabel_output.stdout)

# Test for grep
grep_input = GrepInput(fileInput=receptor, pattern='ATOM')
grep_output = Grep.compute(grep_input)

if debug:
	print("==============================")
	print("GREP OUTPUT:")
	print("==============================")
	print(grep_output.stdout)
