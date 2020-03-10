from qcelemental.models import ProtoModel
from mmelemental.models.util.input import FileInput
from mmelemental.models.chem.codes import ChemCode
from typing import List, Optional
from pydantic import Field

class MMoleculeReaderInput(ProtoModel):
    file: Optional[FileInput] = Field(
        None, 
        description = 'Input filename object.'
    )
    code: Optional[ChemCode] = Field(
        None, 
        description = 'Chemical code object that stores a smiles, smarts, etc. code. See :class: ``Identifiers``.'
    ) 
    resremove: Optional[List[str]] = Field(
        None, 
        description = 'List of residues (by name) to remove.'
    )
    dtype: Optional[str] = Field(
        None,
        description = 'Data type specification e.g. rdkit, gro, pdb, etc. See :class: ``MMoleculeReader``.'
    )