from qcelemental import models
from typing import List, Optional, Union
from pydantic import validator, Field
from pathlib import Path
import os

class CmdOutput(models.ProtoModel):
    stdout_: str = Field(
        ...,
        description = "Standard output."
    )
    stderr_: Optional[str] = Field(
        None,
        description = "Standard error."
    )
    log_: Optional[str] = Field(
        None,
        description = "Logging output"
    )

    class Config(models.ProtoModel.Config):
        fields = {
            "stdout_": "stdout",
            "stderr_": "stderr",
            "log_": "log"
            }

    @property
    def stdout(self):
        return self.stdout_

    @property
    def stderr(self):
        return self.stderr_ 

    @property
    def log(self):
        return self.log_ 

    def dict(self, *args, **kwargs):
        kwargs["by_alias"] = True
        kwargs["exclude_unset"] = True
        return super().dict(*args, **kwargs)

class FileOutput(models.ProtoModel):
    path: str = Field(
        ...,
        description='path to output file.'
    )

    @validator('path')
    def _exists(cls, v):
        if os.path.isfile(v):
            raise IOError(f'File {v} already eixsts.')
        return v

    @property
    def ext(self):
        return Path(self.path).suffix

    def __enter__(self):
        return self

    def write(self, contents: str):
        with open(self.path, 'w') as fp:
            fp.write(contents)

    def remove(self):
        if os.path.isfile(self.path):
            os.remove(self.path)

    def __exit__(self, type, value, tb):
        if not tb:
            self.remove()
        else:
            raise Exception

class ComputeOutput(models.ProtoModel):
    cmdout: Optional[CmdOutput] = Field(
        None,
        description = "Command-line output class which provides stdout, stderr, and log info."
    )