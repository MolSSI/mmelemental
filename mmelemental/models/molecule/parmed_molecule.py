from qcelemental import models
from pydantic import Field, validator
from typing import List, Dict, Any
from .gen_molecule import ToolkitMolecule

try:
    import parmed
except:
    raise ModuleNotFoundError('Make sure rdkit is installed for code validation.')

class ParmedMolecule(ToolkitMolecule):
    mol: parmed.structure.Structure = Field(..., description = 'ParmEd molecule object.')
    dtype: str = Field('parmed', description = 'Data type of mol.')

    @validator('dtype')
    def dtype_static(cls, v):
        assert v == 'parmed', 'dtype for this object must be parmed.'
        return v

    @classmethod
    def build_mol(cls, inputs: Dict[str, Any], dtype: str) -> "ParmedMolecule":
        """ Creates an instance of ParmedMolecule object storing parmed.structure.Structure. 
        This is done by parsing an input file (pdb, gro, ...).
        """
        if inputs.file:
            filename = inputs.file.path
            try:
                pmol = parmed.load_file(filename)
            except:
                raise ValueError(f"File type not supported: {inputs.file.ext}")

        # construct RDKit molecule from identifiers
        elif inputs.code:
            code_type = inputs.code.code_type
            function = getattr(Chem, f"MolFrom{code_type}") # should work since validation already done by ChemCode!
            pmol = function(inputs.code.code)
        else:
            raise ValueError('Missing input file or code.')

        if inputs.code:
            raise NotImplementedError('No support for Chemical codes with ParmEd.')

        return ParmedMolecule(mol=pmol)
