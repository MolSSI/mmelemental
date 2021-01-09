from mmelemental.models.base import Base
from mmelemental.models.util.input import FileInput
from mmelemental.models.util.output import FileOutput
from mmelemental.models.chem.codes import ChemCode
from .gen_molecule import ToolkitMolecule
from typing import Optional, Union, Dict, Any
from pydantic import Field
import importlib

class MolIO(Base):
    _extension_maps = {
        'qcelemental':
        {
            ".npy": "numpy",
            ".json": "json",
            ".xyz": "xyz",
            ".psimol": "psi4",
            ".psi4": "psi4",
            ".msgpack": "msgpack"
        },
        'rdkit':
        {
            ".pdb": "pdb",
            ".mol": "mol",
            ".mol2": "mol2",
            ".tpl": "tpl",
            ".sdf": "sdf",
            ".smi": "smiles"
        },
        'parmed':
        {
            # ".gro": "gro", DOES NOT WORK IN WRITING -> BUG SOMEWHERE?
            ".psf": "psf",
            ".pdb": "pdb",
            ".top": "top"
        }
    }
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
        if args.get('code'):
            if isinstance(args['code'], str):
                args['code'] = ChemCode(code=args['code'])

        super().__init__(**args)

    def files_toolkit(self) -> str:
        """ Returns toolkit name (if any) for reading/writing molecular files. If no
        appropriate toolkit is available on the system, this method raises an error. """

        if hasattr(self, 'file'):
            if self.file:
                dtype_file = self.file.dtype or self.file.ext
            else:
                dtype_file = None
        else:
            dtype_file = None

        if hasattr(self, 'top_file'):
            if self.top_file:
                dtype_top = self.top_file.dtype or self.top_file.ext
            else:
                dtype_top = None
        else:
            dtype_top = None

        if not dtype_file and not dtype_top:
            raise ValueError(f"No files specified to select the appropriate toolkit for in object: {self}.")

        for toolkit in MolInput._extension_maps:
            if dtype_file:
                if MolInput._extension_maps[toolkit].get(dtype_file):
                    if importlib.util.find_spec(toolkit):
                        if dtype_top:
                            if MolInput._extension_maps[toolkit].get(dtype_top):
                                if importlib.util.find_spec(toolkit):
                                    return toolkit
                        else:
                            return toolkit
            else:
                if MolInput._extension_maps[toolkit].get(dtype_top):
                    if importlib.util.find_spec(toolkit):
                        return toolkit

        if dtype_file and dtype_top:
            raise ValueError(f'Could not find appropriate toolkit for reading input files: {self.file.path}, {self.top_file.path}')
        elif dtype_file:
            raise ValueError(f'Could not find appropriate toolkit for reading input file: {self.file.path}')
        else:
            raise ValueError(f'Could not find appropriate toolkit for reading input file: {self.top_file.path}')

class MolInput(MolIO):
    file: Optional[Union[FileInput, str]] = Field(
        None, 
        description = 'Input coords file name or object.'
    ) 
    top_file: Optional[Union[FileInput, str]] = Field(
        None, 
        description = 'Input topology file name or object.'
    )

    def __init__(self, **args):
        file_exists = args.get('file') or args.get('top_file')
        if (file_exists and args.get('code')) or (file_exists and args.get('data')) or (args.get('data') and args.get('code')):
            raise ValueError('Only 1 input type Field (code, file(s), or data) is allowed.')

        if args.get('file'):
            if isinstance(args['file'], str):
                args['file'] = FileInput(path=args['file'])

        if args.get('top_file'):
            if isinstance(args['top_file'], str):
                args['top_file'] = FileInput(path=args['top_file'])

        super().__init__(**args)

class MolOutput(MolIO):
    mol: Any = Field(
        ...,
        description = 'Molecule object such as the :class:``Molecule`` model. '
    )
    file: Optional[Union[FileOutput, str]] = Field(
        None,
        description = 'Output file name.'
    )
    mode: Optional[str] = Field(
        'w',
        description = 'Write mode. By default overwrites existing file. Options: "w" (write) or "a" (append). \
                       See https://docs.python.org/3/tutorial/inputoutput.html#reading-and-writing-files.'
    )

    def __init__(self, **args):
        file_exists = args.get('file')
        if (file_exists and args.get('code')) or (file_exists and args.get('data')) or (args.get('data') and args.get('code')):
            raise ValueError('Only 1 input type Field (code, file, or data) is allowed.')

        if args.get('file'):
            if args.get('mode'):
                mode = args.get('mode')
            else:
                mode = 'w'
            if isinstance(args['file'], str):
                args['file'] = FileOutput(path=args['file'], mode=mode)

        super().__init__(**args)