from mmelemental.models.base import Base
from mmelemental.models.util.input import FileInput
from .gen_ff import ToolkitFF
from typing import Optional, Union, Dict
from pydantic import Field


class FFInput(Base):
    file: Optional[Union[FileInput, str]] = Field(
        None, description="Input params/ff file name or object."
    )
    data: Optional[ToolkitFF] = Field(
        None, description="Toolkit-specific data object e.g. rdkit.Chem.rdchem.Mol"
    )
    args: Optional[Dict] = Field(None, description="Additional arguments to pass.")

    def __init__(self, **args):
        if args.get("file") and args.get("data"):
            raise ValueError("Only 1 input Field (file, or data) is allowed.")

        if args.get("file"):
            if isinstance(args["file"], str):
                args["file"] = FileInput(path=args["file"])

        super().__init__(**args)
