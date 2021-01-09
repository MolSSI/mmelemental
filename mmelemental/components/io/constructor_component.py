from mmelemental.models.base import Base
from mmic.components.blueprints.generic_component import GenericComponent
from typing import Any, Dict, List, Optional, Tuple

from mmelemental.models.molecule.mm_molecule import Molecule
from mmelemental.models.molecule.io_molecule import MolInput
from mmelemental.models.chem.codes import ChemCode
from mmelemental.models.util.input import FileInput

class MolConstructorComponent(GenericComponent):
    """ Class for constructing Molecule from ChemCode or FileInput. """

    @classmethod
    def input(cls):
        return MolInput

    @classmethod
    def output(cls):
        return Molecule

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, Dict[str, Any]]:

        return True, self.constructor(inputs)

    def constructor(self, model: Base) -> Molecule:
        if isinstance(model, ChemCode):
            ctype = str(model.code_type).lower()
            return Molecule(symbols=['C'], geometry=[0,0,0], identifiers={ctype: model})
        elif isinstance(model, FileInput):
            return Molecule.from_file(model.path)
        elif isinstance(model, Molecule):
            return model
        elif isinstance(model, MolInput):
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
        return Molecule

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, Dict[str, Any]]:

        return True, self.constructor(inputs)

    def constructor(self, model: Base) -> Molecule:
        if isinstance(model, ChemCode):
            ctype = str(model.code_type).lower()
            return Molecule(symbols=['C'], geometry=[0,0,0], identifiers={ctype: model})
        elif isinstance(model, FileInput):
            return Molecule.from_file(model.path)
        elif isinstance(model, Molecule):
            return model
        else:
            raise ValueError(f'Input type {type(model)} not supported for {self.__class__}')
