from typing import List, Optional, Any, Dict, Tuple
from mmic.components.blueprints.generic_component import GenericComponent
from mmelemental.models.molecule.io_mol import MolOutput
from mmelemental.components.trans.template_component import TransComponent
from mmelemental.models.base import ToolkitModel
import importlib


class TkMolWriterComponent(GenericComponent):
    @classmethod
    def input(cls):
        return MolOutput

    @classmethod
    def output(cls):
        return ToolkitModel

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, Dict[str, Any]]:

        if isinstance(inputs, dict):
            inputs = TkMolWriterComponent.input()(**inputs)

        if inputs.ext:
            translator = TransComponent.find_molwrite_tk(inputs.ext)

            if not translator:
                raise ValueError(
                    f"Could not write file with ext {inputs.ext}. Please install an appropriate translator."
                )
        else:
            raise ValueError("Data type not supplied.")

        if importlib.util.find_spec(translator):
            mod = importlib.import_module(translator)
            tkmol = mod._classes_map.get("Mol")

        if not tkmol:
            raise ValueError(
                f"No Molecule model found while looking in translator: {translator}."
            )

        return True, tkmol.from_schema(data=inputs.mol)
