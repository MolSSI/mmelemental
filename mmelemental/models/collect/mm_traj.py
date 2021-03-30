from pydantic import Field
from typing import Union, Optional, Tuple, List, Dict, Any
from mmelemental.models.util.input import FileInput
from qcelemental.models.types import Array
from mmelemental.models.molecule.mm_mol import Molecule
from mmelemental.models.base import ProtoModel
from mmelemental.models.util.output import FileOutput
from pathlib import Path
import importlib

from .sm_ensem import Microstate

__all__ = ["Trajectory", "Frame", "TrajInput"]

_trans_nfound_msg = "MMElemental translation requires mmic_translator. \
Solve by: pip install mmic_translator"


class TrajReaderInput(ProtoModel):
    traj: Union[FileInput, str] = Field(..., description="Trajectory input filename.")
    top: Optional[Union[FileInput, str]] = Field(
        ..., description="Topology input filename."
    )


class TrajInput(ProtoModel):
    geometry: Optional[int] = Field(
        ..., description="Atomic positions of length natoms. Default unit is Angstroms."
    )
    geometry_units: Optional[str] = Field(
        "angstrom", description="Units for atomic geometry. Defaults to Angstroms."
    )
    velocities: Optional[int] = Field(
        None,
        description="Atomic velocities of length natoms. Default unit is Angstroms/femtoseconds.",
    )
    velocities_units: Optional[str] = Field(
        "angstrom/fs",
        description="Units for atomic velocities. Defaults to Angstroms/femtoseconds.",
    )
    forces: Optional[int] = Field(
        None, description="Atomic forces of length natoms. KiloJoules/mol.Angstroms."
    )
    forces_units: Optional[str] = Field(
        "kJ/(mol*angstrom)",
        description="Units for atomic forces. Defaults to KiloJoules/mol.Angstroms.",
    )
    freq: Optional[int] = Field(
        None,
        description="Every number of steps the geometry, velocities, and/or forces are sampled.",
    )


class Frame(Microstate):
    timestep: Optional[float] = Field(
        None, description="Timestep size. Default unit is femtoseconds."
    )
    timestep_units: Optional[str] = Field(
        "fs", description="Timestep size units. Defaults to femtoseconds."
    )


class Trajectory(ProtoModel):
    top: Optional[Union[Molecule, List[Molecule]]] = Field(
        None,
        description="Single or multiple :class:``Molecule`` object(s) representing the molecular topology.",
    )
    frames: List[Frame] = Field(
        None, description="A list of :class:``Frame`` objects of length nframes."
    )

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
            in the :class:``TransComponent`` class.
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
            from mmic_translator.components.supported import reg_trans

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
        """Writes the Trajectory to a file.
        Parameters
        ----------
        filename : str
            The filename to write to
        dtype : str, optional
            The type of file to write, attempts to infer dtype from the filename if not provided.
        translator: Optional[str], optional
            Translator name e.g. mmic_rdkit. Takes precedence over dtype. If unset,
            MMElemental attempts to find an appropriate translator if it is registered
            in the :class:``TransComponent`` class.
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

    def to_data(self, dtype: str):
        """ Converts Trajectory to toolkit-specific trajectory object. """
        raise NotImplementedError(f"Data type {dtype} not available.")
