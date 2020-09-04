from mmcomponents.components.blueprints.generic_component import GenericComponent
from mmelemental.models.utils.output import FileOutput

class MoleculeWriter(GenericComponent):

    @classmethod
    def input(cls):
        return Molecule

    @classmethod
    def output(cls):
        return FileOutput