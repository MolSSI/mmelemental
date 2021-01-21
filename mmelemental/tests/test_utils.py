import sys
import os

# import models
from mmelemental.models.util.input import FileInput, OpenBabelInput, GrepInput
from mmelemental.models.chem.codes import ChemCode

# Test input file
receptor = FileInput(path="mmelemental/data/molecules/alanine.pdb")


def test_obabel():
    from mmelemental.components.util.openbabel_component import OpenBabelComponent

    obabel_input = OpenBabelInput(fileInput=receptor, outputExt="pdbqt")
    obabel_output = OpenBabelComponent.compute(obabel_input)


def test_grep():
    from mmelemental.components.util.grep_component import GrepComponent

    grep_input = GrepInput(fileInput=receptor, pattern="ATOM")
    grep_output = GrepComponent.compute(grep_input)
