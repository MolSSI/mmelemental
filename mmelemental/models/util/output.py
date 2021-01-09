from mmelemental.models.base import Base
from typing import List, Optional, Union
from pydantic import validator, Field
from pathlib import Path
import os

class CmdOutput(Base):
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

    class Config(Base.Config):
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

class FileOutput(Base):
    """ Model for writing output to files. No file is created if "write" method is not invoked. """
    path: str = Field(..., description='Output filename path. ')
    clean: bool = Field(False, description='If set to True, the file is removed once object is out of scope.')
    mode: str = Field('w', description='File write mode. Defaults to overwriting files. See https://docs.python.org/3/tutorial/inputoutput.html#reading-and-writing-files.')
    dtype: Optional[str] = Field(None, description='Object data type e.g. PDB, TRR, etc. May not be consistent with the file extension.')

    @validator('path')
    def _exists(cls, v):
        if os.path.isfile(v):
            pass
            # should we: raise IOError(f'File {v} already eixsts.')  ???
        return v

    @property
    def ext(self):
        return Path(self.path).suffix

    @property
    def abs_path(self):
        return os.path.abspath(self.path)

    @property
    def name(self):
        return os.path.basename(self.path)
        
    def __enter__(self):
        return self

    def write(self, contents: str):
        with open(self.abs_path, self.mode) as fp:
            fp.write(contents)

    def remove(self):
        if os.path.isfile(self.abs_path):
            os.remove(self.abs_path)

    @staticmethod
    def rand_name(N=8):
        import string, random
        return ''.join(random.choices(string.ascii_uppercase + string.digits, k=N))

    def __exit__(self, type, value, tb):
        if not tb:
            if self.clean:
                self.remove()
        else:
            raise Exception

class ComputeOutput(Base):
    cmdout: Optional[CmdOutput] = Field(
        None,
        description = "Command-line output class which provides stdout, stderr, and log info."
    )