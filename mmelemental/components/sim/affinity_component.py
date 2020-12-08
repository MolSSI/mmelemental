from mmelemental.models.output.docking import DockingOutput, AffinityOutput
from mmic.components.blueprints.generic_component import GenericComponent

class DockingAffinityComponent(GenericComponent):
    @classmethod
    def input(cls):
        return DockingOutput

    @classmethod
    def output(cls):
        return AffinityOutput
