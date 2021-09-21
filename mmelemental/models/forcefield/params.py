from pydantic import validator, root_validator, Field
from mmelemental.models.base import ProtoModel
from mmelemental.models.util.output import FileOutput
from cmselemental.util.decorators import classproperty
from typing import Any, Optional, Dict, List
import ast
import glob
import pathlib
from importlib.machinery import SourceFileLoader


class Params(ProtoModel):
    """Superclass for different FF sucblasses e.g. Bonds, Angles, etc."""

    _path = None
    form: str = Field(
        ...,
        description="Name or form of the potential. See cls.supported_potentials for available potentials.",
    )
    version: Optional[str] = Field(  # type: ignore
        None,
        description="Version of the force field this model stores. This field can be arbitrary.",
    )
    params: Any = Field(..., description="Specific force field parameters model.")
    defs: Optional[List[str]] = Field(  # type: ignore
        None,
        description="Particle definition. For atomic forcefields, this could be the atom type (e.g. HH31) or SMIRKS (OFF) representation.",
    )
    extras: Dict[str, Any] = Field(  # type: ignore
        None,
        description="Additional information to bundle with the object. Use for schema development and scratch space.",
    )
    # Constructors
    @classmethod
    def from_file(
        cls,
        filename: str,
        dtype: Optional[str] = None,
        translator: Optional[str] = None,
        **kwargs,
    ) -> "Params":
        """
        Constructs a Params object from a file.
        Parameters
        ----------
        filename: str
            The topology or FF filename to build from.
        dtype: Optional[str], optional
            The type of file to interpret e.g. psf. If unset, mmelemental attempts to discover the file type.
        translator: Optional[str], optional
            Translator name e.g. mmic_parmed. Takes precedence over dtype. If unset, MMElemental attempts
            to find an appropriate translator if it is registered in the :class:``TransComponent`` class.
        **kwargs: Optional[Dict[str, Any]], optional
            Any additional keywords to pass to the constructor.
        Returns
        -------
        Params
            A constructed Angles object.
        """
        import json

        fileobj = FileOutput(path=filename)
        dtype = dtype or fileobj.ext.strip(".")

        assert dtype == "json", "Only JSON format supported for now."

        if translator:
            raise NotImplementedError

        with open(filename, "r") as fp:
            data = json.load(fp)

        return cls(**data)

    def to_file(
        self,
        filename: str,
        dtype: Optional[str] = None,
        translator: Optional[str] = None,
        **kwargs: Dict[str, Any],
    ) -> None:
        """Writes the Params to a file.
        Parameters
        ----------
        filename : str
            The filename to write to
        dtype : Optional[str], optional
            The type of file to write (e.g. psf, top, etc.), attempts to infer dtype from
            file extension if not provided.
        translator: Optional[str], optional
            Translator name e.g. mmic_parmed. Takes precedence over dtype. If unset, MMElemental attempts
            to find an appropriate translator if it is registered in the :class:``TransComponent`` class.
        **kwargs: Optional[str, Dict], optional
            Additional kwargs to pass to the constructor.
        """
        if not dtype:
            ext = pathlib.Path(filename).suffix
        else:
            ext = "." + dtype

        if translator or ext != ".json":
            raise NotImplementedError

        stringified = self.json(**kwargs)
        mode = kwargs.pop("mode", "w")

        with open(filename, mode) as fp:
            fp.write(stringified)

    def __eq__(self, other):
        """
        Checks if two molecules are identical. This is a molecular identity defined
        by scientific terms, and not programing terms, so it's less rigorous than
        a programmatic equality or a memory equivalent `is`.
        """

        if isinstance(other, dict):
            other = Params(**other)
        elif isinstance(other, Params):
            pass
        else:
            raise TypeError(f"Comparison not understood of type '{type(other)}'.")

        return self.get_hash() == other.get_hash()

    @validator("form")
    def _registered_name(cls, v):
        if v not in cls.supported_potentials:
            raise NotImplementedError(f"{v} is not supported in MMElemental.")
        return v

    @root_validator
    def _valid_params(cls, values):
        v = values.get("params")
        if isinstance(v, dict):
            v_class = cls.supported_potentials.get(values.get("form"))
            v = v_class(**v)
            values["params"] = v
        assert v.__class__.__name__ == values.get(
            "form"
        ), f"Params type: {v.__class__.__name__} != form: {values.get('form')}."
        return values

    def get_hash(self):
        """
        Returns the hash of the force field object.
        """
        import json, hashlib

        m = hashlib.sha1()
        concat = ""

        # np.set_printoptions(precision=16)
        for field in self.hash_fields:
            data = getattr(self, field)
            if data is not None and field != "params":
                # if field == "nonbonded":
                #    data = float_prep(data, GEOMETRY_NOISE)

                concat += json.dumps(data, default=lambda x: x.ravel().tolist())

        m.update(concat.encode("utf-8"))
        return m.hexdigest()

    # Properties
    @property
    def hash_fields(self):
        return ["params"]

    @property
    def units(self):
        """Returns instance (object) units i.e. any Field name ending with _units."""
        data_units = super().units
        data_units.update(self.params.units)
        return data_units

    @classproperty
    def supported_potentials(cls) -> Dict[str, "Params"]:
        """Returns a dictionary for all available Params subclasses available in cls._path."""
        files = glob.glob(cls._path)
        classes = {}
        for sfile in files:
            with open(sfile, "r") as fp:
                src = fp.read()
                p = ast.parse(src)
                for node in ast.walk(p):
                    if isinstance(node, ast.ClassDef):
                        mod_name = pathlib.Path(sfile).stem
                        mod = SourceFileLoader(mod_name, sfile).load_module()
                        classes[node.name] = getattr(mod, node.name)
        return classes
