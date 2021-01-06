from pydantic import Field, validator
from typing import List, Dict, Any
from .gen_molecule import ToolkitMolecule

try:
    import parmed
except:
    raise ModuleNotFoundError('Make sure parmed is installed for code validation.')

class ParmedMolecule(ToolkitMolecule):
    mol: parmed.structure.Structure = Field(..., description = 'ParmEd molecule object.')

    @property
    def dtype(self):
        return 'parmed'   

    @classmethod
    def build(cls, inputs: Dict[str, Any], dtype: str) -> "ParmedMolecule":
        """
        Creates an instance of ParmedMolecule object storing parmed.structure.Structure. 
        This is done by parsing an input file (pdb, gro, ...).
        """
        if inputs.file:
            coords_fname = inputs.file.path
            if inputs.top_file:
                top_fname = inputs.top_file.path
            else:
                top_fname = None
            try:
                if top_fname:
                    pmol = parmed.load_file(filename=top_fname, xyz=coords_fname)
                else:
                    pmol = parmed.load_file(filename=coords_fname)
            except:
                raise ValueError(f"File type not supported: {inputs.file.ext}")

        elif inputs.code:
            raise NotImplementedError('No support for Chemical codes with ParmEd.')

        return cls(mol=pmol)
