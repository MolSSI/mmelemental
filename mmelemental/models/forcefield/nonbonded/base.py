from mmelemental.models.forcefield.params import Params
import os
import pathlib

__all__ = ["NonBonded"]


class NonBonded(Params):
    """ Model that describes non-bonded interactions between atoms/particles. """

    _path = os.path.join(pathlib.Path(__file__).parent.absolute(), "potentials", "*.py")
