Collections
###########

The :py:mod:`~mmelemental.models.collections` submodule provides models that describe series of states.

Trajectory
==========

.. _traj_descr:

Description
-----------
The `Trajectory` model describes trajectories in classical molecular dynamics. 
Model creation occurs with a kwargs constructor as shown by equivalent operations below:

.. code-block:: python

    >>> traj = mmelemental.models.Trajectory(
            natoms = 1,
            nframes = 3,
            timestep = 1.0,
            timestep_units = "femtoseconds",
            geometry = [2.0, 2.09, 0.0, 2.82, 2.09, 0.58, 1.18, 2.09, 0.58],
            geometry_units = "angstrom",
        )
    >>> traj
     Trajectory(name=None, hash='da39a3e')

In addition, `Trajectory` provides :func:`~mmelemental.models.trajectory.Trajectory.from_data` and :func:`~mmelemental.models.trajectory.Trajectory.from_file` methods to create a trajectory from data (e.g. MDAnalysis.Universe) or file objects. 
The methods :func:`~mmelemental.models.trajectory.Trajectory.to_data` and :func:`~mmelemental.models.trajectory.Trajectory.to_file` enable converting a trajectory to data and file objects, respectively. See the API for more details.
See the :doc:`tutorials_traj` tutorials for more in-depth examples.

.. _trajectory_api:

API
---
.. automodule:: mmelemental.models.collections.trajectory
   :members:
