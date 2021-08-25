Models
######
.. _MMSchema: https://molssi.github.io/mmschema
.. _pydantic: https://sphinx-pydantic.readthedocs.io

MMElemental provides pydantic_ models or data classes based on the MMSchema_ specification. These models
provide serialization and data validation. When suitable, a model hash is automatically generated to allow 
each object to be uniquely identified.


Molecules
=========

The `mmelemental.models.molecule` submodule provides models that describe molecular topology and conformation.

Molecule
--------
The `Molecule` model describes particle systems such as atomistic, coarse-grained, and symbolic (e.g. smiles) molecules. 
Model creation occurs with a kwargs constructor as shown by equivalent operations below:

.. code-block:: python

    >>> mol = mmelemental.models.Molecule(
            **{
                "symbols": ["O", "H", "H"],
                "geometry": [2.0, 2.09, 0.0, 2.82, 2.09, 0.58, 1.18, 2.09, 0.58],
                "connectivity": [(0, 1, 1.0), (0, 2, 1.0)],
            }
        )
    >>> mol
     Molecule(name='H2O', hash='704898f')

In addition, `Molecule` provides ``from_data`` and ``from_file`` methods to create a molecule from data (e.g. MDAnalysis.Universe) or file objects. 
The methods ``to_data`` and ``to_file`` enable converting a molecule to data and file objects, respectively. See the API for more details.


Topology
--------
To be completed soon.

API
---
.. automodule:: mmelemental.models.molecule
   :members:

Force fields
============
The `mmelemental.models.forcefield` submodule provides models that describe force field definitions or parametrized molecules (e.g. data stored in "topology" files).

ForceField
----------
This is the most basic model for storing force field data that describe atomistic or coarse-grained inter-particle potentials.
Model creation occurs with a kwargs constructor as shown by equivalent operations below:

.. code-block:: python

    >>> ff = mmelemental.models.ForceField(
            **{
                "symbols": ["H", "H", "O"],
                "charges": [-0.834, 0.417, 0.417],
                "masses": [16.0, 1.008, 1.008],
                "exclusions": "3",
                "defs": ["OW", "HW", "HW"],
            }
        )
    >>> ff
     ForceField(name='forcefield', form=[], hash='8352e99')

Note that this force field object has no form since we did not define any bonded (or non-bonded) terms as part of the interaction. This can be specified using
additional models explained in the next subsections.

In addition, `ForceField` provides ``from_data`` and ``from_file`` methods to create a forcefield from data (e.g. parmed.Structure) or file objects.
The methods ``to_data`` and ``to_file`` enable converting a force field to data and file objects, respectively. See the API for more details.

NonBonded
---------
To be completed.

Bonds
-----
The `Bonds` model describes pairwise potentials for 2 connected atoms. Model creation occurs with a kwargs constructor as shown by equivalent operations below:

.. code-block:: python

    >>> bonds = mmelemental.models.forcefield.Bonds(
            **{
                "form": "Harmonic",
                "params": {"spring": [2512.08, 2512.08]},
                "lengths": [0.9572, 0.9572],
                "connectivity": [[0,1,1.0], [0,2,1.0]],
            }
        )

This model can be used as an input to `ForceField` for an added descrition of bonded potentials.

.. code-block:: python

    >>> ff = mmelemental.models.ForceField(
            **{
                "symbols": ["H", "H", "O"],
                "charges": [-0.834, 0.417, 0.417],
                "masses": [16.0, 1.008, 1.008],
                "exclusions": "3",
                "defs": ["OW", "HW", "HW"],
                "bonds": bonds,
            }
        )
    >>> ff
     ForceField(name='forcefield', form=['Bonds'], hash='dfbbf4a')

Angles
------
To be completed.

Dihedrals
---------
To be completed.

Improper dihedrals
------------------
To be completed.

API
---
.. automodule:: mmelemental.models.forcefield
   :members:

Collections
===========

Trajectory
----------
To be completed.

Ensemble
--------
To be completed.

API
---
.. automodule:: mmelemental.models.collect
   :members:
