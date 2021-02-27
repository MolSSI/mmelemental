from pydantic import root_validator
from mmelemental.models.base import ProtoModel
import ast
import glob
import os


class Params(ProtoModel):
    _path_name: str

    @root_validator
    def _is_registered(cls, values):
        if cls.__name__ not in cls.supported_potentials():
            raise NotImplementedError(
                f"{cls.__name__} is not supported in MMElemental."
            )
        return values

    def dict(self, *args, **kwargs):
        kwargs["exclude"] = {"provenance"}
        return super().dict(*args, **kwargs)

    @classmethod
    def supported_potentials(cls):
        files = glob.glob(cls._path_name)
        classes = []
        for file in files:
            with open(file, "r") as fp:
                src = fp.read()
                p = ast.parse(src)
                classes.extend(
                    [
                        node.name
                        for node in ast.walk(p)
                        if isinstance(node, ast.ClassDef)
                    ]
                )

        return classes
