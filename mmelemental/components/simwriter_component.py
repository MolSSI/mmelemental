from mmcomponents.components.blueprints.generic_component import GenericComponent
from mmelemental.models.utils.output import FileOutput
from mmelemental.models.sim.base import Base

class SimWriter(GenericComponent):

    @classmethod
    def input(cls):
        return Base

    @classmethod
    def output(cls):
        return FileOutput