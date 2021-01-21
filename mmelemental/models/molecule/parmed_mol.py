from pydantic import Field, validator
from typing import List, Dict, Any
from .gen_mol import ToolkitMol
from mmelemental.util.decorators import require
import parmed

class ParmedMol(ToolkitMol):
    mol: parmed.structure.Structure = Field(..., description = 'ParmEd molecule object.')

    @require('parmed')
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @property
    def dtype(self):
        return 'parmed'   
    
    @classmethod
    @require('parmed')
    def build(cls, inputs: Dict[str, Any], dtype: str) -> "ParmedMol":
        """
        Creates an instance of ParmedMol object storing parmed.structure.Structure. 
        This is done by parsing an input file (pdb, gro, ...).

        .. todo:: use dtype somewhere? Do need it?
        """
        import parmed

        if inputs.file:
            coords_fname = inputs.file.abs_path
            if inputs.top_file:
                top_fname = inputs.top_file.abs_path
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
