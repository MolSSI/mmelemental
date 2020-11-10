from mmcomponents.components.blueprints.generic_component import GenericComponent
from mmelemental.models.util.output import FileOutput
from mmelemental.models.molecule.mda_molecule import MdaMolecule
from mmelemental.models.molecule.mol_reader import MoleculeReaderInput
from mmelemental.models.molecule.mm_molecule import Molecule
from typing import Dict, Any, List, Tuple, Optional

try:
    from MDAnalysis import Universe
except:
    raise ModuleNotFoundError('Make sure MDAnalysis is installed.')

class MoleculeToMda(GenericComponent):
    """ A component for converting Molecule to MDAnalysis molecule object. """
    @classmethod
    def input(cls):
        return Molecule

    @classmethod
    def output(cls):
        return MdaMolecule

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, MdaMolecule]:

        # do stuff here to build mmol from inputs

        return True, MdaMolecule(mol=mmol)

class MdaToMolecule(GenericComponent):
    """ A component for converting MDAnalysis molecule to Molecule object. """

    @classmethod
    def input(cls):
        return MoleculeReaderInput

    @classmethod
    def output(cls):
        return Molecule

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None) -> Tuple[bool, Dict[str, Any]]:
        
        if isinstance(inputs, dict):
            inputs = MoleculeReaderComponent.input()(**inputs)
        
        if inputs.data:
            dtype = inputs.data.dtype
            assert dtype == 'mdanalysis'
            pmol = inputs.data
        elif inputs.code:
            raise NotImplementedError('MDAnalysis does not support instantiating molecule objects from chemical codes.')          
        elif inputs.file:
            dtype = inputs.file.ext
            mmol = MdaMolecule.build_mol(inputs, dtype)
        else:
            raise NotImplementedError(f'Data type {dtype} not yet supported.')
        
        if inputs.args:
            orient = inputs.args.get('orient')
            validate = inputs.args.get('validate')
            kwargs = inputs.args.get('kwargs')
        else:
            orient, validate, kwargs = False, None, None

        # get all properties + more from Universe?

        input_dict = {'symbols': symbs, 
                      'geometry': geo, 
                      'residues': residues, 
                      'connectivity': connectivity,
                      'masses': masses,
                      'names': names}        

        if kwargs:
            input_dict.update(kwargs)

        if inputs.code:
            return True, Molecule(orient=orient, validate=validate, identifiers={dtype: inputs.code}, **input_dict)
        else:
            return True, Molecule(orient=orient, validate=validate, **input_dict)