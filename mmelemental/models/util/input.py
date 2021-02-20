from typing import List, Optional, Union
from mmelemental.models.base import ProtoModel
from pydantic import validator, Field
import os

from .output import FileOutput
from pathlib import Path


class FileInput(ProtoModel):
    """ A model that represents system files that may or may not exist at the time of object insantiation. """

    path: str = Field(..., description="File path, relative or absolute.")
    dtype: Optional[str] = Field(
        None,
        description="Object data type e.g. PDB, TRR, etc. May not be consistent with the file extension.",
    )

    @validator("path")
    def _exists(cls, path):
        if not os.path.isfile(path):
            raise IOError(f"Input file {path} does not exist.")
        return path

    @property
    def abs_path(self):
        return os.path.abspath(self.path)

    @property
    def ext(self):
        return Path(self.path).suffix

    @property
    def name(self):
        return os.path.basename(self.path)

    def __enter__(self):
        return self

    def read(self) -> str:
        with open(self.abs_path, "r") as fp:
            return fp.read()

    def remove(self):
        if os.path.isfile(self.abs_path):
            os.remove(self.abs_path)

    def __exit__(self, type, value, tb):
        if not tb:
            self.remove()
        else:
            raise Exception


class CmdInput(ProtoModel):
    fileInput: Union[FileInput, List[FileInput]] = Field(
        ...,
        description="FileInput object or list of FileInut objects. See the :class: ``FileInput``.",
    )
    fileOutput: Optional[Union[FileOutput, List[FileOutput]]] = Field(
        None,
        description="FileOutput object or list of FileOutput objects. See the :class: ``FileOutput``.",
    )
    args: Optional[List[str]] = Field(
        None, description="List of additional command-line arguments."
    )


class OpenBabelInput(CmdInput):
    outputExt: str = Field(..., description="File output extension.")


class GrepInput(CmdInput):
    pattern: str = Field(..., description="Pattern to search for in input file.")
