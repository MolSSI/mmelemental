from mmelemental.models.input.docking import DockingInput
from mmelemental.models.output.docking import DockingOutput
from mmic.components.blueprints.generic_component import GenericComponent


class DockingComponent(GenericComponent):
    @classmethod
    def input(cls):
        return DockingInput

    @classmethod
    def output(cls):
        return DockingOutput
