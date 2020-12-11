from mmic.components.blueprints.generic_component import GenericComponent
from mmelemental.models.util.output import FileOutput

class MolWriterComponent(GenericComponent):

    @classmethod
    def input(cls):
        return Molecule

    @classmethod
    def output(cls):
        return FileOutput