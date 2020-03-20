from qcelemental.models import ProtoModel
from mmcomponents.components.blueprints.generic_component import GenericComponent
from typing import Any, Dict, List, Optional, Tuple

from mmelemental.models.molecule.mm_molecule import MMolecule
from mmelemental.models.molecule.mol_reader import MMoleculeReaderInput
from mmelemental.models.chem.codes import ChemCode
from mmelemental.models.util.input import FileInput

class MolConstructorComponent(GenericComponent):
    """ Class for constructing MMolecule from ChemCode or FileInput. """

    @classmethod
    def input(cls):
        return MMoleculeReaderInput

    @classmethod
    def output(cls):
        return MMolecule

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, Dict[str, Any]]:

        return True, self.constructor(inputs)

    def constructor(self, model: ProtoModel) -> MMolecule:
        if isinstance(model, ChemCode):
            ctype = str(model.code_type).lower()
            return MMolecule(symbols=['C'], geometry=[0,0,0], identifiers={ctype: model})
        elif isinstance(model, FileInput):
            return MMolecule.from_file(model.path)
        elif isinstance(model, MMolecule):
            return model
        elif isinstance(model, MMoleculeReaderInput):
            if model.code:
                return self.constructor(model.code)
            elif model.file:
                return self.constructor(model.file)
            else:
                raise ValueError('Input file or chemical code must be supplied.')
        else:
            raise ValueError(f'Input type {type(model)} not supported for {self.__class__}')

class ForceFieldConstructorComponent(GenericComponent):
    """ Class for constructing ForceField object from FileInput. """

    @classmethod
    def output(cls):
        return MMolecule

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, Dict[str, Any]]:

        return True, self.constructor(inputs)

    def constructor(self, model: ProtoModel) -> MMolecule:
        if isinstance(model, ChemCode):
            ctype = str(model.code_type).lower()
            return MMolecule(symbols=['C'], geometry=[0,0,0], identifiers={ctype: model})
        elif isinstance(model, FileInput):
            return MMolecule.from_file(model.path)
        elif isinstance(model, MMolecule):
            return model
        else:
            raise ValueError(f'Input type {type(model)} not supported for {self.__class__}')
