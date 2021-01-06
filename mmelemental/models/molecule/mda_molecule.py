from pydantic import Field, validator
from typing import List, Dict, Any
from .gen_molecule import ToolkitMolecule

try:
    import MDAnalysis
except:
    raise ModuleNotFoundError('Make sure MDAnalysis is installed.')

class MdaMolecule(ToolkitMolecule):
    mol: MDAnalysis.Universe = Field(..., description = 'MDAnalysis molecule object.')

    @property
    def dtype(self):
        return 'mdanalysis'   

    @classmethod
    def build(cls, inputs: Dict[str, Any], dtype: str) -> "MdaMolecule":
        """
        Creates an instance of MdaMolecule object storing MDAnalysis.Universe. 
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
                    mmol = MDAnalysis.Universe(top_fname, coords_fname)
                else:
                    mmol = MDAnalysis.Universe(coords_fname)
            except:
                raise ValueError(f"File type not supported: {inputs.file.ext}")

        elif inputs.code:
            raise NotImplementedError('No support for Chemical codes with MDAnalysis.')

        return cls(mol=mmol)
