from mmic.components.blueprints.generic_component import GenericComponent
from typing import Dict, Any, List, Union
import abc
import importlib


class TransComponent(GenericComponent, abc.ABC):
    """ An abstract component that serves as a translator template for converting between MMSchema and other MM codes. """

    _supported_trans = {
        "mmic_mda": "MDAnalysis",
    }
    _supported_versions = {}

    @staticmethod
    def get(obj: object, prop: str) -> Any:
        """ Returns obj.prop if it exists. """
        return getattr(obj, prop) if hasattr(obj, prop) else None

    @staticmethod
    def has(obj: object, prop: str) -> bool:
        """ Returns True if obj.prop exists and is not None. """
        if hasattr(obj, prop):
            return False if getattr(obj, prop) is None else True
        return False

    @staticmethod
    def installed() -> List[str]:
        """ Returns module spec if it exists. """
        return [
            spec
            for spec in TransComponent._supported_trans
            if importlib.util.find_spec(spec)
        ]

    @staticmethod
    def find_molread_ext_maps() -> Dict[str, Dict]:
        """ Returns a Dict of molecule translators and the file formats they can read. """
        trans_mod = (importlib.import_module(mod) for mod in TransComponent.installed())
        return {mod.__name__: mod.models.molread_ext_maps for mod in trans_mod}

    @staticmethod
    def find_molwrite_ext_maps() -> Dict[str, Dict]:
        """ Returns a Dict of molecule translators and the file formats they can write. """
        trans_mod = (importlib.import_module(mod) for mod in TransComponent.installed())
        return {mod.__name__: mod.models.molwrite_ext_maps for mod in trans_mod}

    @staticmethod
    def find_molread_tk(dtype: str) -> Union[str, None]:
        extension_maps = TransComponent.find_molread_ext_maps()
        for toolkit in extension_maps:
            if extension_maps[toolkit].get(dtype):
                if importlib.util.find_spec(toolkit):
                    return toolkit
        return None

    @staticmethod
    def find_molwrite_tk(dtype: str) -> Union[str, None]:
        extension_maps = TransComponent.find_molwrite_ext_maps()
        for toolkit in extension_maps:
            if extension_maps[toolkit].get(dtype):
                if importlib.util.find_spec(toolkit):
                    return toolkit
        return None

    @staticmethod
    def find_trans(dtype: str = None) -> str:
        """Returns mmic_translator name (if any) for writing molecular objects. If no
        appropriate toolkit is available on the system, this method raises an error."""
        for trans, tk in TransComponent._supported_trans.items():
            if dtype == tk:
                return trans

        raise ValueError(f"Could not find appropriate toolkit for {dtype} object.")
