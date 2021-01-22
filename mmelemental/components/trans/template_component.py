from mmic.components.blueprints.generic_component import GenericComponent
from mmelemental.util.decorators import require
from typing import Dict, Any, List, Tuple, Optional


class TransComponent(GenericComponent):
    """ A translator component template for converting between MMSchema and other MM codes. """

    @staticmethod
    def get(obj: object, prop: str) -> Any:
        """ Returns obj.prop if it exists. """
        return getattr(obj, prop) if hasattr(obj, prop) else None

    @staticmethod
    def has(obj: object, prop: str) -> bool:
        """ Returns True if obj.prop exists and is not None. """
        if hasattr(obj, prop):
            return False if getattr(obj, prop) is None else True
        else:
            False
