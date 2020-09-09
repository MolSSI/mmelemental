from mmelemental.models.base import Base
from pydantic import Field
from typing import Any, Array

class MolecularDynamics(Base):
	""" Molecular dynamics parameter input schema."""

	# System fields
	mol: "MMolecule" = Field(..., description = 'Molecule mechanics object.')
	
	# I/O fields
	traj: str = Field(None, description = "Trajectory filename for coordinates and/or velocities, forces.")
	traj_vel: str = Field(None, description = "Trajectory filename explicitly for velocities.")
	traj_for: str = Field(None, description = "Trajectory filename explicitly for forces.")
	traj_freq: int = Field(None, description = "Frequency of writing to trajectory output file.")
	restart: str = Field(None, description = "Restart filename.")
	restart_freq: int = Field(None, description = "Frequency of writing to restart output file.")

	# Thermodynamic fields
    temp: Array[float] = Field(None, description = 'Temperature(s) for thermostat(s) in Kelvin.')
    press: Array[float] = Field(None, description = 'Pressure(s) for barostate(s) in Bar.')
