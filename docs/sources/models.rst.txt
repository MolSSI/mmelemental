Models
######
.. _MMSchema: https://molssi.github.io/mmschema
.. _pydantic: https://sphinx-pydantic.readthedocs.io

MMElemental provides models based on the MMSchema_ specification. These models
provide serialization and data validation based on the pydantic_ library. Each model
has a set of fields that, when suitable, are used to automatically generate a unique hash 
that enables each object to be uniquely identified.


Molecules
=========

The :py:mod:`~mmelemental.models.molecule` submodule provides models that describe molecular topology and conformation.

Topology
--------
To be completed soon.


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

In addition, `Molecule` provides :func:`~mmelemental.models.molecule.Molecule.from_data` and :func:`~mmelemental.models.molecule.Molecule.from_file` methods to create a molecule from data (e.g. MDAnalysis.Universe) or file objects. 
The methods :func:`~mmelemental.models.molecule.Molecule.to_data` and :func:`~mmelemental.models.molecule.Molecule.to_file` enable converting a molecule to data and file objects, respectively. See the API for more details.


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
                "defs": ["OW", "HW", "HW"],
            }
        )
    >>> ff
     ForceField(name='forcefield', form=[], hash='8352e99')

Note that this force field object has no form since we did not define any bonded (or non-bonded) terms as part of the interaction. This can be specified using
additional models explained in the next subsections.

In addition, `ForceField` provides :func:`~mmelemental.models.forcefield.ForceField.from_data` and :func:`~mmelemental.models.forcefield.ForceField.from_file` methods to create a forcefield from data (e.g. parmed.Structure) or file objects.
The methods :func:`~mmelemental.models.forcefield.ForceField.to_data` and :func:`~mmelemental.models.forcefield.ForceField.to_file` enable converting a force field to data and file objects, respectively. See the API for more details.

NonBonded
---------
The `NonBonded` model describes potentials for non-connected particles. Model creation occurs with a kwargs constructor as shown by equivalent operations below:

.. code-block:: python

    >>> nb = mmelemental.models.forcefield.NonBonded(
            **{
                "form": "LennardJones",
                "params": {"epsilon": [0.636386, 0.0, 0.0], "sigma": [1.575305, 0.0, 0.0]},
            }
        )

This model can be used as an input to `ForceField` for an extended description of non-bonded interactions.

.. code-block:: python

    >>> ff = mmelemental.models.ForceField(
            **{
                "symbols": ["H", "H", "O"],
                "charges": [-0.834, 0.417, 0.417],
                "masses": [16.0, 1.008, 1.008],
                "defs": ["OW", "HW", "HW"],
                "nonbonded": nb,
            }
        )
    >>> ff
     ForceField(name='forcefield', form=['NonBonded'], hash='0cbf0de')

Bonds
-----
The `Bonds` model describes pairwise potentials for 2 connected particles. Model creation occurs with a kwargs constructor as shown by equivalent operations below:

.. code-block:: python

    >>> bonds = mmelemental.models.forcefield.Bonds(
            **{
                "form": "Harmonic",
                "params": {"spring": [2512.08, 2512.08]},
                "lengths": [0.9572, 0.9572],
                "connectivity": [[0,1,1.0], [0,2,1.0]],
            }
        )

This model can be used as an input to `ForceField` for an extended description of chemical bonds.

.. code-block:: python

    >>> ff = mmelemental.models.ForceField(
            **{
                "symbols": ["H", "H", "O"],
                "charges": [-0.834, 0.417, 0.417],
                "masses": [16.0, 1.008, 1.008],
                "exclusions": "3",
                "defs": ["OW", "HW", "HW"],
                "nonbonded": nb,
                "bonds": bonds,
            }
        )
    >>> ff
     ForceField(name='forcefield', form=['NonBonded', 'Bonds'], hash='236e292')

Angles
------
The `Angles` model describes 3-body angle potentials for 3 connected particles. Model creation occurs with a kwargs constructor as shown by equivalent operations below:

.. code-block:: python

    >>> angles = mmelemental.models.forcefield.Angles(
            **{
                "form": "Harmonic",
                "params": {"spring": [0.09565291598722435]},
                "angles": [104.52],
                "connectivity": [
                    [0, 1, 2],
                ],
            }
        )

This model can be used as an input to `ForceField` for an extended description of chemical bonds.

.. code-block:: python

    >>> ff = mmelemental.models.ForceField(
            **{
                "symbols": ["H", "H", "O"],
                "charges": [-0.834, 0.417, 0.417],
                "masses": [16.0, 1.008, 1.008],
                "exclusions": "3",
                "defs": ["OW", "HW", "HW"],
                "nonbonded": nb,
                "bonds": bonds,
                "angles": angles,
            }
        )
    >>> ff
     ForceField(name='forcefield', form=['NonBonded', 'Bonds', 'Angles'], hash='6bdf578')

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
