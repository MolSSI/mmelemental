from pydantic import Field
from typing import Union, Optional, Tuple, Dict, Any
from mmelemental.models.util.input import FileInput
from qcelemental.models.types import Array
from mmelemental.models.molecule.mm_molecule import Molecule
from mmelemental.models.base import Base
from mmelemental.components.trajreader_component import SingleFrameComponent, MultiFrameComponent

class TrajectoryReaderInput(Base):
    traj: Union[FileInput, str] = Field(
        ...,
        description = 'Trajectory input filename.'
    )
    top: Optional[Union[FileInput, str]] = Field(
        ...,
        description = 'Topology input filename.'
    )

class TrajectoryOutput(Base):
    coordinates: Array[float] = Field(
        ..., 
        description = 'Atomic coordinates of length natoms.'
    )
    velocities: Optional[Array[float]] = Field(
        None, 
        description = 'Atomic velocities of length natoms.'
    )
    forces: Optional[Array[float]] = Field(
        None, 
        description = 'Atomic forces of length natoms.'
    )
    pot_energy: Optional[Array[float]] = Field(
        None, 
        description = 'Total system potential energy.'
    )
    timestep: Optional[float] = Field(
        None, 
        description = 'Timestep size.'
    )
    score: Optional[float] = Field(
        None, 
        description = 'Ligand score used in docking simulations.'
    )

class Trajectory(Base):
    topology: Optional[Array[Molecule]] = Field(
        None, 
        description = 'One or more Molecule object(s).'
    )
    result: Array[TrajectoryOutput] = Field(
        None,
        description = 'An Array of TrajectoryOutput objects of length nframes.'
    )
    _formats: Dict[str, Tuple[str]] =  {
            'dcd': ('mdanalysis', 'mdtraj', 'pytraj', 'loos'),
            'netcdf3': ('mdanalysis', 'mdtraj', 'pytraj', 'loos', 'parmed'),
            'netcdf4': ('mdanalysis', 'mdtraj', 'pytraj', 'loos', 'parmed'),
            'trr': ('mdanalysis', 'mdtraj', 'pytraj', 'loos'),
            'xtc': ('mdanalysis', 'mdtraj', 'pytraj', 'loos')
    }

    @property
    def formats(self) -> Dict[str, Tuple[str]]:
        return self._formats

    # Constructors
    @classmethod
    def from_file(
        cls, 
        traj: Union[FileInput, str],
        top: Union[FileInput, str] = None, 
        dtype: str = None,
        *, 
        all_frames: bool = False,
        **kwargs
    ) -> "Trajectory":
        """
        Constructs a Trajectory object from an input file.
        Parameters
        ----------
        traj: FileInput or str
            Trajectory file to construct object from
        top: FileInput or str, optional
            Topology file to read
        dtype : str, optional
            The type of file to interpret. If not set, mmelemental attempts to discover the file type.
        all_frames: bool, optional
            Reads all frames at once.
        **kwargs: Dict[str, Any]
            Additional kwargs to pass to the constructors. kwargs take precedence over data.
        Returns
        -------
        Trajectory
            A constructed Trajectory class.
        """
        traj_input = TrajectoryReaderInput(traj=traj, top=top)

        if all_frames:
            return MultiFrameComponent.compute(traj_input)
        else:
            return SingleFrameComponent.compute(traj_input)

    @classmethod
    def from_data(cls, 
        data: Any, 
        dtype: Optional[str] = None,
        **kwargs: Dict[str, Any]
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
        if not dtype:
            if data.__class__.__name__ == "Universe":
                from mmelemental.components.mda_component import UniverseToTrajectory
                return UniverseToTrajectory.compute(data)
            else:
                raise NotImplementedError(f'Data type {dtype} not supported by mmelemental.')
        elif dtype == 'mdanalysis':
            from mmelemental.components.mda_component import UniverseToTrajectory
            return UniverseToTrajectory.compute(data)
        else:
            raise NotImplementedError(f'Data type {dtype} not supported by mmelemental.')

    def to_file(self, filename: str, dtype: Optional[str] = None) -> None:
        """ Writes the Trajectory to a file.
        Parameters
        ----------
        filename : str
            The filename to write to
        dtype : str, optional
            The type of file to write, attempts to infer dtype from the filename if not provided.
        """
        raise NotImplementedError(f'Data type {dtype} not available.')

    def to_data(self, dtype: str):
        """ Converts Trajectory to toolkit-specific trajectory object. """
        raise NotImplementedError(f'Data type {dtype} not available.')
