from qcengine.util import execute
from mmelemental.models.util.output import FileOutput
from mmic.components.blueprints import SpecificComponent
from typing import Any, Dict, List, Tuple, Optional, Union
from mmelemental.models.util.output import CmdOutput
from mmelemental.models.util.input import CmdInput


class CmdComponent(SpecificComponent):
    """ Cmd process: build_input() -> run() -> parse_output() -> [clean()] """

    @classmethod
    def input(cls):
        return CmdInput

    @classmethod
    def output(cls):
        return CmdOutput

    def clean(self, files: Union[List[FileOutput], FileOutput]):
        if isinstance(files, list):
            for file in files:
                file.remove()
        else:
            files.remove()

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, Dict[str, Any]]:

        execute_input = self.build_input(inputs)
        exe_success, proc = self.run(execute_input)

        if exe_success:
            out = True, self.parse_output(proc, inputs)
            if execute_input.get("clean_files"):
                self.clean(execute_input.get("clean_files"))
            return out
        else:
            raise ValueError(proc["stderr"])

    def run(
        self,
        inputs: Dict[str, Any],
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, Dict[str, Any]]:

        infiles = inputs["infiles"]
        outfiles = inputs["outfiles"]

        if extra_outfiles is not None:
            outfiles.extend(extra_outfiles)

        command = inputs["command"]
        if extra_commands is not None:
            command.extend(extra_commands)

        exe_success, proc = execute(
            command,
            infiles=infiles,
            outfiles=outfiles,
            scratch_directory=inputs["scratch_directory"],
            scratch_name=scratch_name,
            timeout=timeout,
            environment=inputs.get("environment", None),
        )

        return exe_success, proc

    def build_input(
        self,
        input_model: Dict[str, Any],
        config: Optional["TaskConfig"] = None,
        template: Optional[str] = None,
    ) -> Dict[str, Any]:
        raise NotImplementedError

    def parse_output(
        self, output: Dict[str, str], inputs: Dict[str, Any]
    ) -> Dict[str, Any]:
        raise NotImplementedError
