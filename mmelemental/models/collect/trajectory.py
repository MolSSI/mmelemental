from pydantic import Field, validator, root_validator, constr
from typing import Union, Optional, List, Dict, Any
from cmselemental.types import Array
from mmelemental.models.struct.topology import Topology
from mmelemental.models.base import ProtoModel, Provenance, provenance_stamp
from mmelemental.models.util.output import FileOutput
from pathlib import Path
import importlib
import hashlib
import numpy
from mmelemental.util.units import TIME_DIM, MASS_DIM, LENGTH_DIM, SUBS_DIM
from mmelemental.util.data import (
    float_prep,
    NUMPY_INT,
    NUMPY_FLOAT,
    GEOMETRY_NOISE,
    VELOCITY_NOISE,
    FORCE_NOISE,
)
from cmselemental.util import which_import
import json


__all__ = ["Trajectory"]

_trans_nfound_msg = "MMElemental translation requires mmic_translator. \
Solve by: pip install mmic_translator"


mmschema_trajectory_default = "mmschema_trajectory"


class Trajectory(ProtoModel):
    """Representation of trajectories in classical mechanics. By default, this model loads a single frame in memory."""

    # Basic field
    schema_name: constr(
        strip_whitespace=True, regex=mmschema_trajectory_default
    ) = Field(  # type: ignore
        mmschema_trajectory_default,
        description=(
            f"The MMSchema specification to which this model conforms. Explicitly fixed as {mmschema_trajectory_default}."
        ),
    )
    schema_version: Optional[int] = Field(  # type: ignore
        1,
        description="The version number of ``schema_name`` to which this model conforms.",
    )
    name: Optional[str] = Field(  # type: ignore
        None,
        description="Common or human-readable name to assign to this molecule. This field can be arbitrary; see "
        "``identifiers`` for well-defined labels.",
    )
    # Definition fields
    timestep: Union[Array[float], float] = Field(
        ..., description="Timestep size. Default unit is femtoseconds."
    )
    timestep_units: Optional[str] = Field(
        "fs",
        description="Timestep size units. Defaults to femtoseconds.",
        dimensionality=TIME_DIM,
    )
    natoms: Union[int, Array[numpy.dtype(NUMPY_INT)]] = Field(
        ...,
        description="Number of atoms. A scalar for closed systems and an array for open (variable N) systems.",
    )
    nframes: int = Field(..., description="Number of frames.")  # type: ignore

    ndim: Optional[int] = Field(  # type: ignore
        3, description="Number of spatial dimensions."
    )

    # For storing topological data or time-dependent topologies (e.g. reactive ffs)
    top: Optional[Union[Topology, List[Topology]]] = Field(  # type: ignore
        None,
        description=Topology.__doc__.split("\n\n\n").pop(0),
    )
    # Particle dynamical fields
    geometry: Optional[Array[numpy.dtype(NUMPY_FLOAT)]] = Field(  # type: ignore
        None,
        description="An ordered (natom*ndim*nframes,) array for XYZ atomic coordinates. Default unit is Angstrom. "
        "Storage is sequential in each dimension:\n"
        "[\n"
        "x1o,...,xno, y1o,...,yno, z1o,...,zno,\n"
        "...,\n"
        "x1f,...,xnf, y1f,...,ynf, z1f,...,znf,\n"
        "]",
    )
    geometry_units: Optional[str] = Field(  # type: ignore
        "angstrom",
        description="Units for atomic geometry. Defaults to Angstroms.",
        dimensionality=LENGTH_DIM,
    )
    velocities: Optional[Array[numpy.dtype(NUMPY_FLOAT)]] = Field(  # type: ignore
        None,
        description="An ordered (natoms*ndim*nframes,) array for XYZ atomic velocities. Default unit is "
        "Angstroms/femtoseconds."
        "Storage is sequential in each dimension:\n"
        "[\n"
        "x1o,...,xno, y1o,...,yno, z1o,...,zno,\n"
        "...,\n"
        "x1f,...,xnf, y1f,...,ynf, z1f,...,znf,\n"
        "]",
    )
    velocities_units: Optional[str] = Field(  # type: ignore
        "angstrom/fs",
        description="Units for atomic velocities. Defaults to Angstroms/femtoseconds.",
        dimensionality=LENGTH_DIM / TIME_DIM,
    )
    forces: Optional[Array[numpy.dtype(NUMPY_FLOAT)]] = Field(  # type: ignore
        None,
        description="An ordered (natoms*ndim*nframes,) array for XYZ atomic forces. Default unit is "
        "kJ/mol*angstrom."
        "Storage is sequential in each dimension:\n"
        "[\n"
        "x1o,...,xno, y1o,...,yno, z1o,...,zno,\n"
        "...,\n"
        "x1f,...,xnf, y1f,...,ynf, z1f,...,znf,\n"
        "]",
    )
    forces_units: Optional[str] = Field(  # type: ignore
        "kJ/(mol*angstrom)",
        description="Units for atomic forces. Defaults to kJ/mol*angstrom.",
        dimensionality=MASS_DIM * LENGTH_DIM / (SUBS_DIM * TIME_DIM ** 2),
    )
    provenance: Provenance = Field(
        provenance_stamp(__name__),
        description="The provenance information about how this object (and its attributes) were generated, "
        "provided, and manipulated.",
    )
    extras: Dict[str, Any] = Field(  # type: ignore
        None,
        description="Additional information to bundle with this object. Use for schema development and scratch space.",
    )

    class Config(ProtoModel.Config):
        repr_style = lambda self: [("name", self.name), ("hash", self.get_hash()[:7])]
        fields = {
            "nframes_": "nframes",
        }

    def __repr_args__(self) -> "ReprArgs":
        return [("name", self.name), ("hash", self.get_hash()[:7])]

    def __hash__(self):
        return hash(self.get_hash())

    @validator("geometry", "velocities", "forces")
    def _must_be_n(cls, v, values):
        if v is None:
            return v

        natoms = values["natoms"]
        nframes = values["nframes"]
        try:
            v.reshape(natoms, values["ndim"], nframes)
        except (ValueError, AttributeError):
            raise ValueError("Array must be castable to shape (natoms,ndim,frames)!")

        if v.ndim > 1:
            assert v.ndim == 2, f"Array must be either 1D or 2D! Given ndim: {v.ndim}."
            assert (
                v.shape[1] == values["ndim"] * nframes
            ), f"Array column number ({v.shape[1]}) should be equal to ndim * nframes ({values['ndim'] * nframes})."

        assert values["ndim"] * natoms * nframes == len(
            v
        ), "Number of frames is inconsistent with array shape!"
        return v

    # Properties
    @property
    def hash_fields(self):
        return [
            "geometry",
            "geometry_units",
            "velocities",
            "velocities_units",
            "forces",
            "forces_units",
            "top",
            "timestep",
            "timestep_units",
            "nframes",
            "ndim",
        ]

    def get_geometry(self, frame: int):
        """Returns geometry at a specific snapshot/frame.

        Parameters
        ----------
        frame: int
            Frame number ranges from 0 ... nframes-1

        Returns
        -------
        Array[float]
            Geometry 1D numpy array of length natoms * ndim

        """
        if frame >= self.nframes:
            raise ValueError("Frame number cannot exceed total number of frames!")

        if self.geometry is not None:
            nfree = self.ndim * self.natoms
            return self.geometry[frame * nfree : (frame + 1) * nfree]
        return None

    # Helper methods
    def get_hash(self):
        """Returns the hash of a single frame in a trajectory."""

        m = hashlib.sha1()
        concat = ""

        for field in self.hash_fields:
            data = getattr(self, field)
            if data is not None:
                if field == "geometry":
                    data = float_prep(data, GEOMETRY_NOISE)
                elif field == "velocities":
                    data = float_prep(data, VELOCITY_NOISE)
                elif field == "forces":
                    data = float_prep(data, FORCE_NOISE)

                concat += json.dumps(
                    data, default=lambda x: x.ravel().tolist()
                )  # if serialization fails, assume type is numpy.ndarray

        m.update(concat.encode("utf-8"))
        return m.hexdigest()

    # Constructors
    @classmethod
    def from_file(
        cls,
        traj_filename: str,
        top_filename: Optional[str] = None,
        dtype: str = None,
        *,
        translator: Optional[str] = None,
        all_frames: bool = False,
        **kwargs,
    ) -> "Trajectory":
        """
        Constructs a Trajectory object from an input file.

        Parameters
        ----------
        traj_filename : str
            The atomic positions filename to read
        top_filename: str, optional
            The topology i.e. connectivity filename to read
        dtype : str, optional
            The type of file to interpret. If not set, mmelemental attempts to discover the file type.
        translator: Optional[str], optional
            Translator name e.g. mmic_rdkit. Takes precedence over dtype. If unset,
            MMElemental attempts to find an appropriate translator if it is registered
            in the :class:`TransComponent` class.
        all_frames: bool, optional
            Reads all frames in memory.
        **kwargs: Dict[str, Any], optional
            Additional kwargs to pass to the constructors.

        Returns
        -------
        Trajectory
            A constructed Trajectory class.

        """
        file_ext = Path(traj_filename).suffix
        dtype = dtype or file_ext.strip(".")

        if file_ext in [".json"]:
            if not all_frames:
                raise ValueError(
                    f"Single frame cannot be read from {file_ext} file. Use all_frames=True instead."
                )

            if top_filename:
                raise TypeError(
                    "Molecule topology must be supplied in a single JSON (or similar) file."
                )

            import json

            # Raw string type, read and pass through
            if dtype == "json":
                with open(traj_filename, "r") as infile:
                    data = json.load(infile)
                dtype = "dict"
            else:
                raise KeyError(f"Data type not supported: {dtype}.")

            return cls.from_data(data, dtype=dtype, **kwargs)

        fileobj = FileOutput(path=traj_filename)
        top_fileobj = FileOutput(path=top_filename) if top_filename else None

        # Generic translator component
        try:
            from mmic_translator.components import TransComponent
        except Exception:
            TransComponent = None

        if not translator:
            if not TransComponent:
                raise ModuleNotFoundError(_trans_nfound_msg)
            from mmic_translator import reg_trans

            reg_trans = list(reg_trans)

            while not translator:
                translator = TransComponent.find_trajread_tk(file_ext, trans=reg_trans)
                if not translator:
                    raise ValueError(
                        f"Could not read traj file with ext {file_ext}. Please install an appropriate translator."
                    )
                # Make sure we can import the translator module
                if importlib.util.find_spec(translator):
                    mod = importlib.import_module(translator)

                # If top if supplied, make sure the translator supports the top file extension
                if top_fileobj:
                    top_ext = top_fileobj.ext
                    if top_ext not in mod.ffread_ext_maps:
                        reg_trans.remove(translator)
                        translator = None
                    if not len(reg_trans):
                        raise ValueError(
                            f"Could not read traj and top files with exts {file_ext} and {top_ext}. \
                            Please install an appropriate translator."
                        )
        elif importlib.util.find_spec(translator):
            mod = importlib.import_module(translator)

        tktraj_class = mod._classes_map.get("Trajectory")

        if not tktraj_class:
            raise ValueError(
                f"No Trajectory model found while looking in translator: {translator}."
            )

        tk_traj = tktraj_class.from_file(
            filename=fileobj.abs_path if fileobj else None,
            top_filename=top_fileobj.abs_path if top_fileobj else None,
            dtype=dtype,
            all_frames=all_frames,
        )

        return cls.from_data(tk_traj, dtype=tk_traj.dtype, **kwargs)

    @classmethod
    def from_data(
        cls, data: Any, dtype: Optional[str] = None, **kwargs: Dict[str, Any]
    ) -> "Trajectory":
        """
        Constructs a Trajectory object from a data object.

        Parameters
        ----------
        data: Any
            Data to construct Molecule from
        dtype: str, optional
            How to interpret the data. If not set, mmelemental attempts to discover this based on input type.
            Possible values: mdanalysis, mdraj, pytraj
        **kwargs: Dict[str, Any]
            Additional kwargs to pass to the constructors. kwargs take precedence over data.

        Returns
        -------
        Trajectory
            A constructed Trajectory class.

        """
        if isinstance(data, dict):
            kwargs.pop("dtype", None)  # remove dtype if supplied
            kwargs.update(data)
            return cls(**kwargs)

        return data.to_schema(**kwargs)

    def to_file(
        self,
        filename: str,
        dtype: Optional[str] = None,
        *,
        translator: Optional[str] = None,
        **kwargs: Optional[Dict[str, Any]],
    ) -> None:
        """
        Writes the Trajectory to a file.

        Parameters
        ----------
        filename : str
            The filename to write to
        dtype : str, optional
            The type of file to write, attempts to infer dtype from the filename if not provided.
        translator: Optional[str], optional
            Translator name e.g. mmic_rdkit. Takes precedence over dtype. If unset,
            MMElemental attempts to find an appropriate translator if it is registered
            in the :class:`TransComponent` class.
        **kwargs: Optional[Dict[str, Any]], optional
            Additional kwargs to pass to the constructor.

        """
        if not dtype:
            from pathlib import Path

            ext = Path(filename).suffix
        else:
            ext = "." + dtype

        mode = kwargs.pop("mode", "w")

        if ext == ".json":
            stringified = self.json(**kwargs)
            with open(filename, mode) as fp:
                fp.write(stringified)
        elif ext == ".yaml" or ext == ".yml":
            stringified = self.yaml(**kwargs)
            with open(filename, mode) as fp:
                fp.write(stringified)
        else:  # look for an installed mmic_translator
            try:
                from mmic_translator.components import TransComponent
            except Exception:
                TransComponent = None

            if not TransComponent:
                raise ModuleNotFoundError(_trans_nfound_msg)
            translator = TransComponent.find_molwrite_tk(ext)

            if not translator:
                raise NotImplementedError(
                    f"File extension {ext} not supported with any installed translators."
                )

            tk_traj = self.to_data(translator=translator, **kwargs)
            tk_traj.to_file(filename, dtype=dtype, **kwargs)  # pass dtype?

    def to_data(
        self,
        dtype: Optional[str] = None,
        *,
        translator: Optional[str] = None,
        **kwargs: Optional[Dict[str, Any]],
    ) -> "ToolkitModel":
        """
        Converts Molecule to toolkit-specific molecule (e.g. rdkit, MDAnalysis, parmed).

        Parameters
        ----------
        dtype: str, optional
            The type of data object to convert to e.g. mdanalysis, rdkit, parmed, etc.
        translator: Optional[str], optional
            Translator name e.g. mmic_rdkit. Takes precedence over dtype. If unset,
            MMElemental attempts to find an appropriate translator if it is registered
            in the :class:`TransComponent` class.
        **kwargs: Optional[Dict[str, Any]], optional
            Additional kwargs to pass to the constructor.

        Returns
        -------
        ToolkitModel
            Toolkit-specific molecule model

        """
        try:
            from mmic_translator.components import TransComponent
        except Exception:
            TransComponent = None

        if not translator:
            if not TransComponent:
                raise ModuleNotFoundError(_trans_nfound_msg)
            if not dtype:
                raise ValueError(
                    f"Either translator or dtype must be supplied when calling {__name__}."
                )
            translator = TransComponent.find_trans(dtype)

        if which_import(translator, return_bool=True):
            mod = importlib.import_module(translator)
            tkmol = mod._classes_map.get("Trajectory")

            if not tkmol:
                raise ValueError(
                    f"No Molecule model found while looking in translator: {translator}."
                )

            return tkmol.from_schema(self)
        else:
            raise NotImplementedError(
                f"translator {translator} not available. Make sure it is properly installed."
            )
