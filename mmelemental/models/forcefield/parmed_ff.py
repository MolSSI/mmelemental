from pydantic import Field, validator
from typing import List, Dict, Any
from .gen_ff import ToolkitFF
from .io_ff import FFInput

try:
    import parmed
except:
    raise ModuleNotFoundError("Make sure parmed is installed for code validation.")


class ParmedFF(ToolkitFF):
    ff: parmed.structure.Structure = Field(
        ..., description="ParmEd force field object."
    )

    @property
    def dtype(self):
        return "parmed"

    @classmethod
    def build(cls, inputs: FFInput, dtype: str) -> "ParmedFF":
        """
        Creates an instance of ParmedFF object storing a subclass of parmed.structure.Structure.
        E.g. parmed.charmm.psf.CharmmPsfFile. This is done by parsing an input file (psf, top, ...).

        .. todo:: guess dtype?
        """
        if not isinstance(inputs, FFInput):
            inputs = FFInput(**inputs)

        if inputs.file:
            top_fname = inputs.file.abs_path
            ff = parmed.load_file(filename=top_fname)
        else:
            raise NotImplementedError(
                "No support for anything but input files for now."
            )
        return cls(ff=ff)
