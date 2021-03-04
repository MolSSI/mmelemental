from mmelemental.models.forcefield.params import Params
import os
import pathlib

__all__ = ["DihedralsParams"]


class DihedralsParams(Params):
    _path_name = os.path.join(
        pathlib.Path(__file__).parent.absolute(), "potentials", "*.py"
    )
