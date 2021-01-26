from mmelemental.components.util.cmd_component import CmdComponent
from typing import Any, Dict, Optional
import os
from mmelemental.models.util.input import OpenBabelInput
from mmelemental.models.util.output import CmdOutput


class OpenBabelComponent(CmdComponent):
    @classmethod
    def input(cls):
        return OpenBabelInput

    @classmethod
    def output(cls):
        return CmdOutput

    def build_input(
        self,
        input_model: OpenBabelInput,
        config: "TaskConfig" = None,
        template: Optional[str] = None,
    ) -> Dict[str, Any]:
        output = "tmp." + input_model.outputExt
        cmd = ["obabel", input_model.fileInput.abs_path, "-O" + output]

        if input_model.args:
            for arg in input_model.args:
                cmd.append(arg)

        env = os.environ.copy()

        if config:
            env["MKL_NUM_THREADS"] = str(config.ncores)
            env["OMP_NUM_THREADS"] = str(config.ncores)

        scratch_directory = config.scratch_directory if config else None

        return {
            "command": cmd,
            "infiles": None,
            "outfiles": [output],
            "scratch_directory": scratch_directory,
            "environment": env,
        }

    def parse_output(
        self, outfiles: Dict[str, Dict[str, str]], input_model: OpenBabelInput
    ) -> CmdOutput:
        output_file = outfiles["outfiles"]["tmp." + input_model.outputExt]

        return CmdOutput(stdout=output_file)
