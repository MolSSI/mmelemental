from mmelemental.components.util.cmd_component import CmdComponent
from typing import Any, Dict, List, Optional, Tuple
import os
from models.components.utils.input import GrepInput
from models.components.utils.output import CmdOutput

class Grep(CmdComponent):

    @classmethod
    def input(cls):
        return GrepInput

    @classmethod
    def output(cls):
        return CmdOutput

    def build_input(
        self, inputs: Dict[str, Any], config: Optional["TaskConfig"] = None, template: Optional[str] = None
    ) -> Dict[str, Any]:
        
        cmd = ["grep"]

        args = inputs.args
        input_model = {'input': inputs.fileInput.path, 'pattern': inputs.pattern, 'args': args}

        if input_model['args']:
            for arg in input_model['args']:
                cmd.append(arg)

        cmd.append(input_model['pattern']) 

        if isinstance(input_model['input'], list):
            for ginput in input_model['input']:
                cmd.append(ginput)
        elif isinstance(input_model['input'], str):
            cmd.append(input_model['input'])
        else:
            raise Exception

        env = os.environ.copy()

        if config:
            env["MKL_NUM_THREADS"] = str(config.ncores)
            env["OMP_NUM_THREADS"] = str(config.ncores)

        scratch_directory = config.scratch_directory if config else None

        return {
            "command": cmd,
            "infiles": None,
            "outfiles": None,
            "scratch_directory": scratch_directory,
            "environment": env
        }

    def parse_output(self, outfiles: Dict[str, str], input_model: GrepInput) -> CmdOutput:
        
        output_file = outfiles['stdout']

        return CmdOutput(stdout=output_file)