from typing import List, Optional, Any, Dict, Tuple

from mmic.components.blueprints.generic_component import GenericComponent
from mmelemental.models.molecule.mol_reader import MolReaderInput
from mmelemental.models.molecule.gen_molecule import ToolkitMolecule
import qcelemental

class MolReaderComponent(GenericComponent):
    """ Factory component that constructs a Molecule object from MolReaderInput.
    Which toolkit-specific component is called depends on MolReaderInput.data.dtype."""

    @classmethod
    def input(cls):
        return MolReaderInput

    @classmethod
    def output(cls):
        from mmelemental.models.molecule.mm_molecule import Molecule
        return Molecule

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None) -> Tuple[bool, Dict[str, Any]]:
        
        from mmelemental.models.molecule.mm_molecule import Molecule

        if isinstance(inputs, dict):
            inputs = MolReaderComponent.input()(**inputs)

        if inputs.args:
            orient = inputs.args.get('orient')
            validate = inputs.args.get('validate')
            kwargs = inputs.args.get('kwargs')
        else:
            orient, validate, kwargs = False, None, None

        if inputs.data:
            dtype = inputs.data.dtype
            if dtype == 'qcelem':
                qmol = qcelemental.models.molecule.Molecule.from_data(data, dtype, orient=orient, validate=validate, **kwargs)
                return True, Molecule(orient=orient, validate=validate, **qmol.to_dict())
            elif dtype == 'rdkit':
                from mmelemental.components.trans.rdkit_component import RDKitToMolecule
                return True, RDKitToMolecule.compute(inputs)
            elif dtype == 'parmed':
                from mmelemental.components.trans.parmed_component import ParmedToMolecule
                return True, ParmedToMolecule.compute(inputs)               
            else:
                raise NotImplementedError(f'Data type not yet supported: {dtype}.')
        # Only RDKit is handling chem codes and file objects for now!
        elif inputs.code:
            from mmelemental.components.trans.rdkit_component import RDKitToMolecule
            return True, RDKitToMolecule.compute(inputs)
        elif inputs.file:
            from mmelemental.components.trans.rdkit_component import RDKitToMolecule
            return True, RDKitToMolecule.compute(inputs)
        else:
            raise NotImplementedError('Molecules can be instantiated from codes, files, or other data objects.')

class TkMolReaderComponent(GenericComponent):

    _extension_maps = {
        'qcelem':
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
            ".smiles": "smiles"
        },
        'parmed':
        {
            ".gro": "gro",
            ".psf": "psf",
            ".pdb": "pdb",
            ".top": "top"
        }
    }

    @classmethod
    def input(cls):
        return MolReaderInput

    @classmethod
    def output(cls):
        return ToolkitMolecule

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None) -> Tuple[bool, Dict[str, Any]]:
        
        if isinstance(inputs, dict):
            inputs = TkMolReaderComponent.input()(**inputs)

        if inputs.file:
            for ext_map_key in TkMolReaderComponent._extension_maps:
                dtype = TkMolReaderComponent._extension_maps[ext_map_key].get(inputs.file.ext)
                if dtype:
                    if inputs.top_file:
                        if TkMolReaderComponent._extension_maps[ext_map_key].get(inputs.top_file.ext):
                            break
                    else:
                        break

        elif inputs.code:
            dtype = inputs.code.code_type.lower()
        else:
            # need to support TkMolecule conversion from e.g. rdkit to parmed, etc.
            raise ValueError('Data type not understood. Supply a file or a chemical code.')

        if '.' + dtype in TkMolReaderComponent._extension_maps['rdkit']:
            from mmelemental.models.molecule.rdkit_molecule import RDKitMolecule
            return True, RDKitMolecule.build(inputs, dtype)
        elif '.' + dtype in TkMolReaderComponent._extension_maps['parmed']:
            from mmelemental.models.molecule.parmed_molecule import ParmedMolecule
            return True, ParmedMolecule.build(inputs, dtype)
        else:
            raise ValueError(f'Data type {dtype} not supported by {self.__class__}.')
