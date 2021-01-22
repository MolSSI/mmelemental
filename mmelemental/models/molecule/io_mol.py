from mmelemental.models.base import Base
from mmelemental.models.util.input import FileInput
from mmelemental.models.util.output import FileOutput
from mmelemental.models.chem.codes import ChemCode
from .gen_mol import ToolkitMol
from typing import Optional, Union, Dict, List, Any
from pydantic import Field
import importlib
import abc


class Translators:
    """ A wrapper class that provides methods for discovering available/installed MMSchema translators. """

    @staticmethod
    def supported() -> Dict[str, str]:
        return {
            "mmic_mda": "MDAnalysis",
        }

    @staticmethod
    def installed() -> List[str]:
        """ Returns module spec if it exists. """
        return [
            spec for spec in Translators.supported() if importlib.util.find_spec(spec)
        ]

    @staticmethod
    def find_molread_ext_maps() -> Dict[str, Dict]:
        """ Returns a Dict of molecule translators and the file formats they can read. """
        trans_mod = (importlib.import_module(mod) for mod in Translators.installed())
        return {mod.__name__: mod.models.molread_ext_maps for mod in trans_mod}

    @staticmethod
    def find_molwrite_ext_maps() -> Dict[str, Dict]:
        """ Returns a Dict of molecule translators and the file formats they can write. """
        trans_mod = (importlib.import_module(mod) for mod in Translators.installed())
        return {mod.__name__: mod.models.molwrite_ext_maps for mod in trans_mod}

    @staticmethod
    def find_molread_tk(dtype: str) -> Union[str, None]:
        extension_maps = Translators.find_molread_ext_maps()
        for toolkit in extension_maps:
            if extension_maps[toolkit].get(dtype):
                if importlib.util.find_spec(toolkit):
                    return toolkit
        return None

    @staticmethod
    def find_molwrite_tk(dtype: str) -> Union[str, None]:
        extension_maps = Translators.find_molwrite_ext_maps()
        for toolkit in extension_maps:
            if extension_maps[toolkit].get(dtype):
                if importlib.util.find_spec(toolkit):
                    return toolkit
        return None

    @staticmethod
    def find_trans(dtype: str = None) -> str:
        """Returns mmic_translator name (if any) for writing molecular objects. If no
        appropriate toolkit is available on the system, this method raises an error."""
        for trans, tk in Translators.supported().items():
            if dtype == tk:
                return trans

        raise ValueError(f"Could not find appropriate toolkit for {dtype} object.")


class MolIO(Base, abc.ABC):
    code: Optional[Union[ChemCode, str]] = Field(
        None,
        description="Chemical code object that stores a smiles, smarts, etc. code. See :class:``Identifiers``.",
    )
    data: Optional[ToolkitMol] = Field(
        None, description="Toolkit-specific data object e.g. :class:``MdaMolecule``."
    )
    kwargs: Optional[Dict] = Field(None, description="Additional arguments to pass.")

    def __init__(self, **args):
        if args.get("code"):
            if isinstance(args["code"], str):
                args["code"] = ChemCode(code=args["code"])

        super().__init__(**args)

    @abc.abstractmethod
    def find_toolkit(self) -> str:
        """Must return toolkit name (if any) for reading/writing molecular files. If no
        appropriate toolkit is available on the system, this method must raise an error.
        See Implementation in :class:``MolInput`` and :class:``MolOutput``."""
        raise NotImplementedError


class MolInput(MolIO):
    """A model used as an input for constructing a :class:``Mol`` object from file or other data objects.
    See for example :class:``TkMolReaderComponent``."""

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
            ".smi": "smiles",
        },
        "parmed": {
            # ".gro": "gro", DOES NOT WORK IN WRITING -> BUG SOMEWHERE?
            ".psf": "psf",
            ".pdb": "pdb",
            ".top": "top",
        },
        "MDAnalysis": {
            ".gro": "gro",
            ".pdb": "pdb",
            ".top": "top",
            ".psf": "psf",
        },
    }

    file: Optional[Union[FileInput, str]] = Field(
        None, description="Input coords file name or object."
    )
    top_file: Optional[Union[FileInput, str]] = Field(
        None, description="Input topology file name or object."
    )

    def __init__(self, **args):
        file_exists = args.get("file") or args.get("top_file")
        if (
            (file_exists and args.get("code"))
            or (file_exists and args.get("data"))
            or (args.get("data") and args.get("code"))
        ):
            raise ValueError(
                "Only 1 input type Field (code, file(s), or data) is allowed."
            )

        if args.get("file"):
            if isinstance(args["file"], str):
                args["file"] = FileInput(path=args["file"])

        if args.get("top_file"):
            if isinstance(args["top_file"], str):
                args["top_file"] = FileInput(path=args["top_file"])

        super().__init__(**args)

    def find_toolkit(self) -> str:
        """Returns toolkit name (if any) for reading molecular files. If no
        appropriate toolkit is available on the system, this method raises an error."""
        tk_file, tk_top, dtype_file, dtype_top = None, None, None, None

        if self.file:
            dtype_file = self.file.dtype or self.file.ext
            tk_file = Translators.find_toolkit(dtype_file)

        if self.top_file:
            dtype_top = self.top_file.dtype or self.top_file.ext
            tk_top = Translators.find_toolkit(dtype_top)

        if not tk_file and not tk_top:
            raise ValueError(
                f"Could not find appropriate toolkit for reading input files: {self.file.path}, {self.top_file.path}"
            )
        elif not tk_file and not dtype_top:
            raise ValueError(
                f"Could not find appropriate toolkit for reading input file: {self.file.path}"
            )
        elif not tk_top and not dtype_file:
            raise ValueError(
                f"Could not find appropriate toolkit for reading input file: {self.top_file.path}"
            )

        if tk_top and tk_file:
            return tk_top
        else:
            return tk_top or tk_file


class MolOutput(MolIO):
    """A model used as input for converting :class:``Mol`` to other data objects.
    See for example :class:``TkMolWriterComponent``."""

    mol: Any = Field(
        ..., description="Input molecule object such as the :class:``Mol`` model. "
    )
    dtype: Optional[str] = Field(None, description="Output data type.")

    def __init__(self, **args):
        file_exists = args.get("file")
        if (
            (file_exists and args.get("code"))
            or (file_exists and args.get("data"))
            or (args.get("data") and args.get("code"))
        ):
            raise ValueError(
                "Only 1 input type Field (code, file, or data) is allowed."
            )

        if args.get("file"):
            if args.get("mode"):
                mode = args.get("mode")
            else:
                mode = "w"
            if isinstance(args["file"], str):
                args["file"] = FileOutput(path=args["file"], mode=mode)

        super().__init__(**args)

    def find_toolkit(self, dtype: str = None) -> str:
        """Returns toolkit name (if any) for writing molecular files. If no
        appropriate toolkit is available on the system, this method raises an error."""
        if not dtype:
            dtype = self.dtype

        tk_file = Translators.find_molwrite_tk(dtype)

        if not tk_file:
            raise ValueError(
                f"Could not find appropriate toolkit for writing a {dtype} object."
            )

        return tk_file
