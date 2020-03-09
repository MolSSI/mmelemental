from typing import List, Optional, Union
from qcelemental import models
from pydantic import validator, Field
import os

from .output import FileOutput
from pathlib import Path

class FileInput(models.ProtoModel):
    path: str

    @validator('path')
    def _exists(cls, v):
        if not os.path.isfile(v):
            raise IOError(f'Input file {v} does not eixst.')

        return v

    @property
    def ext(self):
        return Path(self.path).suffix

    def read(self) -> str:
        with open(self.path, 'r') as fp:
            return fp.read()

class CmdInput(models.ProtoModel):
    fileInput: Union[FileInput, List[FileInput]] = Field(
        ..., 
        description = 'FileInput object or list of FileInut objects. See the :class: ``FileInput``.'
    )
    fileOutput: Optional[Union[FileOutput, List[FileOutput]]] = Field(
        None,
        description = 'FileOutput object or list of FileOutput objects. See the :class: ``FileOutput``.'
    )
    args: Optional[List[str]] = Field(
        None,
        description = 'List of additional command-line arguments.'
    )

class OpenBabelInput(CmdInput):
    outputExt: str = Field(
        ..., 
        description = 'File output extension.'
    )

class GrepInput(CmdInput):
    pattern: str = Field(
        ...,
        description = 'Pattern to search for in input file.'
    )