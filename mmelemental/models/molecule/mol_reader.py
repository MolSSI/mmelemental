from mmelemental.models.base import Base
from mmelemental.models.util.input import FileInput
from mmelemental.models.chem.codes import ChemCode
from mmelemental.models.molecule.gen_molecule import ToolkitMolecule
from typing import Optional, Union, Dict
from pydantic import Field

class MolReaderInput(Base):
    file: Optional[Union[FileInput, str]] = Field(
        None, 
        description = 'Input coords file name or object.'
    ) 
    top_file: Optional[Union[FileInput, str]] = Field(
        None, 
        description = 'Input topology file name or object.'
    )    
    code: Optional[Union[ChemCode, str]] = Field(
        None, 
        description = 'Chemical code object that stores a smiles, smarts, etc. code. See :class: ``Identifiers``.'
    )
    data: Optional[ToolkitMolecule] = Field(
        None,
        description = 'Toolkit-specific data object e.g. rdkit.Chem.rdchem.Mol'
    ) 
    args: Optional[Dict] = Field(
        None,
        description = 'Additional arguments to pass qcelemental molecule.'
    )

    def __init__(self, **args):
        if (args.get('file') and args.get('code')) or (args.get('file') and args.get('data')) or (args.get('data') and args.get('code')):
            raise ValueError('Only 1 input Field (code, file, or data) is allowed.')

        if args.get('file'):
            if isinstance(args['file'], str):
                args['file'] = FileInput(path=args['file'])

        if args.get('code'):
            if isinstance(args['code'], str):
                args['code'] = ChemCode(code=args['code'])

        super().__init__(**args)