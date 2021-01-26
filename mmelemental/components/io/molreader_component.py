from typing import List, Optional, Any, Dict, Tuple
from mmic.components.blueprints.generic_component import GenericComponent
from mmelemental.models.molecule.io_mol import MolInput
from mmelemental.components.trans.template_component import TransComponent
from mmelemental.models.base import ToolkitModel
import qcelemental
import importlib


class TkMolReaderComponent(GenericComponent):

    _extension_maps = {
        "qcelemental": {
            ".npy": "numpy",
            ".json": "json",
            ".xyz": "xyz",
            ".psimol": "psi4",
            ".psi4": "psi4",
            ".msgpack": "msgpack",
        },
        "rdkit": {
            ".pdb": "pdb",
            ".mol": "mol",
            ".mol2": "mol2",
            ".tpl": "tpl",
            ".sdf": "sdf",
            ".smiles": "smiles",
        },
        "parmed": {
            ".gro": "gro",
            ".psf": "psf",
            ".pdb": "pdb",
            ".top": "top",
            ".sdf": "sdf",
            ".mol": "mol",
            ".mol2": "mol2",
        },
        "MDAnalysis": {
            ".gro": "gro",
            ".pdb": "pdb",
            ".top": "top",
            ".psf": "psf",
        },
    }

    @classmethod
    def input(cls):
        return MolInput

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
            inputs = TkMolReaderComponent.input()(**inputs)

        if inputs.file:
            translator = TransComponent.find_molread_tk(inputs.file.ext)

            if not translator:
                raise ValueError(
                    f"Could not read file with ext {inputs.file.ext}. Please install an appropriate translator."
                )
        else:
            raise ValueError("Data type not understood. Supply a file.")

        if importlib.util.find_spec(translator):
            mod = importlib.import_module(translator + ".models")
            tkmol = mod.classes_map.get("Mol")

        if not tkmol:
            raise ValueError(
                f"No Molecule model found while looking in translator: {translator}."
            )

        return True, tkmol.from_file(
            filename=inputs.file.abs_path if inputs.file else None,
            top_filename=inputs.top_file.abs_path if inputs.top_file else None,
            dtype=inputs.dtype,
        )
