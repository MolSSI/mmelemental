from mmcomponents.components.blueprints.generic_component import GenericComponent
from mmelemental.models.util.output import FileOutput
from mmelemental.models.sim.base import Base
from typing import Dict, List, Any, Optional, Tuple

class SimWriter(GenericComponent):

    @classmethod
    def input(cls):
        return Base

    @classmethod
    def output(cls):
        return FileOutput

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, Dict[str, Any]]:

    	filename = inputs.filename

    	with open(filename, 'w') as fp:
    		fp.write('TESTING')

    	return True, FileOutput(path=filename)