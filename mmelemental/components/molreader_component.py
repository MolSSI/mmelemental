try:
    from rdkit import rdBase, Chem
    from rdkit.Chem import AllChem
except:
    raise ModuleNotFoundError('Make sure rdkit is properly installed on your system.')

from qcelemental.models import ProtoModel
from typing import List, Optional, Any, Dict, Tuple

from mmcomponents.components.blueprints.generic_component import GenericComponent
from mmelemental.models.molecule.molreader import MMoleculeReaderInput
from mmelemental.models.molecule.gen_molecule import ToolkitMolecule


class MMoleculeReader(GenericComponent):

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
            ".code": "code"
        },
        'parmed':
        {
            ".gro": "gro"
        }
    }

    @classmethod
    def input(cls):
        return MMoleculeReaderInput

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
            inputs = MMoleculeReader.input()(**inputs)

        dtype = inputs.dtype

        if not dtype: # try to guess data type
            if inputs.file:
                for ext_map_key in MMoleculeReader._extension_maps:
                    dtype = MMoleculeReader._extension_maps[ext_map_key].get(inputs.file.ext)
                    if dtype:
                        break

            elif inputs.code:
                dtype = '.code'
            else:
                raise ValueError('Data type not understood. Supply a file or a chemical code.')

        if '.' + dtype in MMoleculeReader._extension_maps['rdkit']:
            from mmelemental.models.molecule.rdkit_molecule import RDKitMolecule
            return True, RDKitMolecule.build_mol(inputs, dtype)
        elif '.' + dtype in MMoleculeReader._extension_maps['parmed']:
            from mmelemental.models.molecule.parmed_molecule import ParmedMolecule
            return True, ParmedMolecule.build_mol(inputs, dtype)
        else:
            raise ValueError(f'Data type {dtype} not supported by {self.__class__}.')