Molecules
#########

The :py:mod:`~mmelemental.models.molecule` submodule provides models that describe molecular topology and conformation.

Topology
========
To be completed soon.


Molecule
========

.. _molecule_descr:

Description
-----------
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
See the :doc:`tutorials_mol` tutorials for more in-depth examples.

.. _molecule_api:

API
---
.. automodule:: mmelemental.models.molecule
   :members:
