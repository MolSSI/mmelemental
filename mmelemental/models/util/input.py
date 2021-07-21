from typing import List, Optional, Dict
from mmelemental.models.base import ProtoModel
from pydantic import validator, Field
import os
from pathlib import Path

__all__ = ["FileInput", "CmdInput"]


class FileInput(ProtoModel):
    """A model that represents system files that may or may not exist at the time of object insantiation."""

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
    file_input: List[str] = Field(
        ...,
        description="Input file path(s) or name(s).",
    )
    file_output: Optional[List[str]] = Field(
        None,
        description="Output file path(s) or name(s).",
    )
    flags: Optional[List[str]] = Field(None, description="List of command-line flags.")
    kwargs: Optional[Dict[str, str]] = Field(
        None, description="List of additional command-line arguments."
    )


class OpenBabelInput(CmdInput):
    outputExt: str = Field(..., description="File output extension.")
