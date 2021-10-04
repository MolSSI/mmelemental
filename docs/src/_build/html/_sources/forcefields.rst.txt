Force fields
############
The `mmelemental.models.forcefield` submodule provides models that describe force field definitions or parametrized molecules (e.g. data stored in "topology" files).

ForceField
==========

.. _ffdescr:

Description
-----------
This is the most basic model for storing force field data that describe atomistic or coarse-grained inter-particle potentials.
Model creation occurs with a kwargs constructor as shown by equivalent operations below:

.. code-block:: python

    >>> ff = mmelemental.models.ForceField(
            **{
                "name": "water_ff",
                "symbols": ["O", "H", "H"],
                "charges": [-0.834, 0.417, 0.417],
                "masses": [16.0, 1.008, 1.008],
                "defs": ["OW", "HW", "HW"],
            }
        )
    >>> ff
     ForceField(name='water_ff', form=[], hash='8352e99')

Note that this force field object has no form since we did not define any bonded (or non-bonded) terms as part of the interaction. This can be specified using
additional models explained in the next subsections.

In addition, `ForceField` provides :func:`~mmelemental.models.forcefield.ForceField.from_data` and :func:`~mmelemental.models.forcefield.ForceField.from_file` methods to create a forcefield from data (e.g. parmed.Structure) or file objects.
The methods :func:`~mmelemental.models.forcefield.ForceField.to_data` and :func:`~mmelemental.models.forcefield.ForceField.to_file` enable converting a force field to data and file objects, respectively. See the API for more details.

.. _ffapi:

API
---
.. automodule:: mmelemental.models.forcefield.forcefield
   :members:


NonBonded
=========

.. _nbdescr:

Description
-----------
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
                "name": "water_ff",
                "symbols": ["O", "H", "H"],
                "charges": [-0.834, 0.417, 0.417],
                "masses": [16.0, 1.008, 1.008],
                "defs": ["OW", "HW", "HW"],
                "nonbonded": nb,
            }
        )
    >>> ff
     ForceField(name='water_ff', form=['NonBonded'], hash='0cbf0de')

.. _nbapi:

API
---
.. automodule:: mmelemental.models.forcefield.nonbonded
   :members:

Bonds
=====

.. _bondsdescr:

Description
-----------
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
                "name": "water_ff",
                "symbols": ["O", "H", "H"],
                "charges": [-0.834, 0.417, 0.417],
                "masses": [16.0, 1.008, 1.008],
                "exclusions": "3",
                "defs": ["OW", "HW", "HW"],
                "nonbonded": nb,
                "bonds": bonds,
            }
        )
    >>> ff
     ForceField(name='water_ff', form=['NonBonded', 'Bonds'], hash='236e292')

.. _bondsapi:

API
---
.. automodule:: mmelemental.models.forcefield.bonded.bonds
   :members:

Angles
======

.. _anglesdescr:

Description
-----------
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

This model can be used as an input to `ForceField` for an extended description of 3-angle terms in the force field.

.. code-block:: python

    >>> ff = mmelemental.models.ForceField(
            **{
                "name": "water_ff",
                "symbols": ["O", "H", "H"],
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
     ForceField(name='water_ff', form=['NonBonded', 'Bonds', 'Angles'], hash='6bdf578')

.. _anglesapi:

API
---
.. automodule:: mmelemental.models.forcefield.bonded.angles
   :members:

Dihedrals
---------
To be completed.

Improper dihedrals
------------------
To be completed.

.. _dihedralsapi:

API
---
.. automodule:: mmelemental.models.forcefield.bonded.dihedrals
   :members:
