from mmelemental.models.base import Base
from mmcomponents.components.blueprints.generic_component import GenericComponent
from typing import Any, Dict, List, Optional, Tuple

from mmelemental.models.output.trajectory import Trajectory, TrajectoryReaderInput

class SingleFrameComponent(GenericComponent):
    """ Class for constructing Trajectory that reads a single frame from file(s). """

    @classmethod
    def input(cls):
        return TrajectoryReaderInput

    @classmethod
    def output(cls):
        return Trajectory

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, Dict[str, Any]]:

        return True, Trajectory(topology=..., result=...)

class MultiFrameComponent(GenericComponent):
    """ Class for constructing Trajectory that reads all frames from file(s). """

    @classmethod
    def input(cls):
        return TrajectoryReaderInput

    @classmethod
    def output(cls):
        return Trajectory

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, Dict[str, Any]]:

        return True, Trajectory(topology=..., result=...)