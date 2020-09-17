try:
    from rdkit import rdBase, Chem
    from rdkit.Chem import AllChem
except:
    raise ModuleNotFoundError('Make sure rdkit is properly installed on your system.')

from qcelemental.models import ProtoModel
from typing import List, Optional, Any, Dict, Tuple

from mmcomponents.components.blueprints.generic_component import GenericComponent
from mmelemental.models.molecule.mol_reader import MoleculeReaderInput
from mmelemental.models.molecule.gen_molecule import ToolkitMolecule

class MoleculeReaderComponent(GenericComponent):

    @classmethod
    def input(cls):
        return MoleculeReaderInput

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
        
        # we need to add: mmelemental.models.molecule.parmed_molecule import ParmedMolecule
        from mmelemental.models.molecule.mm_molecule import Molecule

        if isinstance(inputs, dict):
            inputs = MoleculeReaderComponent.input()(**inputs)

        if inputs.args:
            orient = inputs.args.get('orient')
            validate = inputs.args.get('validate')
            kwargs = inputs.args.get('kwargs')
        else:
            orient, validate, kwargs = False, None, None

        if inputs.data:
            dtype = inputs.data.obj_type

            if dtype == 'rdkit':
                rdmol = inputs.data
            else:
                # convert inputs.data to Molecule
                qmol = qcelemental.models.molecule.Molecule.from_data(data, dtype, orient=orient, validate=validate, **kwargs)
                return True, Molecule(orient=orient, validate=validate, **qmol.to_dict())
        elif inputs.code:
            dtype = inputs.code.code_type.lower()
            from mmelemental.models.molecule.rdkit_molecule import RDKitMolecule
            rdmol = RDKitMolecule.build_mol(inputs, dtype)            
        elif inputs.file:
            dtype = inputs.file.ext
            from mmelemental.models.molecule.rdkit_molecule import RDKitMolecule
            rdmol = RDKitMolecule.build_mol(inputs, dtype)
        else:
            raise NotImplementedError
        
        from mmelemental.models.molecule.rdkit_molecule import Bond

        symbs = [atom.GetSymbol() for atom in rdmol.mol.GetAtoms()]

        try:
            residues = [(atom.GetPDBResidueInfo().GetResidueName(), 
                        atom.GetPDBResidueInfo().GetResidueNumber()) 
                        for atom in rdmol.mol.GetAtoms()]
            names = [atom.GetPDBResidueInfo().GetName() for atom in rdmol.mol.GetAtoms()]
        except:
            residues = None
            names = None
        
        connectivity = []

        for bond in rdmol.mol.GetBonds():
            bondOrder = Bond.orders.index(bond.GetBondType())
            connectivity.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), 1)) # replace 1 with bondOrder, 
            # for now this is a hack for qcelemental does not allow bond orders higher than 5

        geo = rdmol.mol.GetConformer(0).GetPositions()

        input_dict = {'symbols': symbs, 
                      'geometry': geo, 
                      'residues': residues, 
                      'connectivity': connectivity,
                      'names': names}        

        if kwargs:
            input_dict.update(kwargs)

        if inputs.code:
            return True, Molecule(orient=orient, validate=validate, identifiers={dtype: inputs.code}, **input_dict)
        else:
            return True, Molecule(orient=orient, validate=validate, **input_dict)

class TkMoleculeReaderComponent(GenericComponent):

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
            ".gro": "gro"
        }
    }

    @classmethod
    def input(cls):
        return MoleculeReaderInput

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
            inputs = TkMoleculeReaderComponent.input()(**inputs)

        if inputs.file:
            for ext_map_key in TkMoleculeReaderComponent._extension_maps:
                dtype = TkMoleculeReaderComponent._extension_maps[ext_map_key].get(inputs.file.ext)
                if dtype:
                    break

        elif inputs.code:
            dtype = inputs.code.code_type.lower()
        else:
            # need to support TkMolecule conversion from e.g. rdkit to parmed, etc.
            raise ValueError('Data type not understood. Supply a file or a chemical code.')

        if '.' + dtype in TkMoleculeReaderComponent._extension_maps['rdkit']:
            from mmelemental.models.molecule.rdkit_molecule import RDKitMolecule
            return True, RDKitMolecule.build_mol(inputs, dtype)
        elif '.' + dtype in TkMoleculeReaderComponent._extension_maps['parmed']:
            from mmelemental.models.molecule.parmed_molecule import ParmedMolecule
            return True, ParmedMolecule.build_mol(inputs, dtype)
        else:
            raise ValueError(f'Data type {dtype} not supported by {self.__class__}.')
