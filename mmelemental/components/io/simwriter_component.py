from mmic.components.blueprints.generic_component import GenericComponent
from mmelemental.models.util.output import FileOutput
from mmelemental.models.input.sim_writer import SimWriterInput
from typing import Dict, List, Any, Optional, Tuple

class SimWriter(GenericComponent):

    @classmethod
    def input(cls):
        return SimWriterInput

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
    ) -> Tuple[bool, FileOutput]:

        file = FileOutput(path=inputs.filename)
        schema = inputs.model
        for key, value in schema:

            if key == 'mol':
                value.to_file(file.name + '.pdb')
                value.to_file(file.name + '.psf')

                if 'NAMD' in inputs.engine:
                    file.write('coordinates \t {}\n'.format(file.name + '.pdb'))
                    file.write('structure \t {}\n'.format(file.name + '.psf'))

            elif value and key != 'cell' and key != 'forcefield' and key != 'solvent' and key != 'provenance':
                file.write('{} \t {}\n'.format(key, value))

        return True, file
