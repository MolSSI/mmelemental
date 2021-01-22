from typing import List, Optional, Any, Dict, Tuple
from mmic.components.blueprints.generic_component import GenericComponent
from mmelemental.models.molecule.io_mol import MolInput
from mmelemental.components.trans.template_component import TransComponent
from mmelemental.models.molecule.gen_mol import ToolkitMol
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
        return ToolkitMol

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
            toolkit = TransComponent.find_molread_tk(inputs.file.ext)

            if not toolkit:
                raise ValueError(
                    f"Data type not understood for file ext {inputs.file.ext}."
                )

        elif inputs.code:
            dtype = inputs.code.code_type.lower()
            toolkit = "rdkit"  # need to support more toolkits for handling chem codes
        else:
            raise ValueError(
                "Data type not understood. Supply a file or a chemical code."
            )

        if toolkit == "rdkit":
            from mmelemental.models.molecule.rdkit_mol import RDKitMol

            return True, RDKitMol.build(inputs, dtype)
        elif toolkit == "parmed":
            from mmelemental.models.molecule.parmed_mol import ParmedMol

            return True, ParmedMol.build(inputs, dtype)
        elif toolkit == "mmic_mda":
            from mmic_mda.models import MdaMol

            if inputs.top_file and inputs.file:
                return True, MdaMol.from_file(
                    filename=inputs.file.abs_path,
                    top_filename=inputs.top_file.abs_path,
                    dtype=inputs.dtype,
                )
            elif inputs.top_file:
                return True, MdaMol.from_file(
                    top_filename=inputs.top_file.abs_path, dtype=inputs.dtype
                )
            elif inputs.file:
                return True, MdaMol.from_file(
                    filename=inputs.file.abs_path, dtype=inputs.dtype
                )
            else:
                raise TypeError("No file was supplied!")
        else:
            raise ValueError(f"Data type {dtype} not supported by {self.__class__}.")
