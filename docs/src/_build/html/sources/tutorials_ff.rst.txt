Force fields
############
.. _MMSchema: https://molssi.github.io/mmschema
.. _pydantic: https://sphinx-pydantic.readthedocs.io
.. _ForceField: forcefields.html#forcefield
.. _Molecule: molecules.html#molecule
.. _Topology: molecules.html#topology
.. _Jupyter: https://jupyter.org
.. _nglview: http://nglviewer.org/nglview/latest/
.. _pint: https://pint.readthedocs.io
.. _MMIC: https://github.com/MolSSI/mmic

The ForceField_ model in MMElemental represents the functional form of force fields and their parameter sets (or subsets for parameterized molecules). ForceField_ objects
can be used in molecular mechanics, molecular dynamics, Monte Carlo methods, molecular docking, coarse-graining, and any other method that utilizes force fields.


Instantiation 
-------------

ForceField_ objects can be created by supplying keyword arguments to the constructor. Per the MMSchema_ specification, there are no required fields. Hence, the simplest way to create a force field object is via:

.. code-block:: python

    >>> ff = mmelemental.models.ForceField()

This creates a Molecule_ object with a unique hash code.

.. code-block:: python

    >>> ff
     ForceField(name='forcefield', form=[], hash='96b52d9')

We can also assign a unique name to this force field object:

.. code-block:: python

    >>> ff = mmelemental.models.ForceField(name="my_forcefield")
    >>> ff
     ForceField(name='my_forcefield', form=[], hash='96b52d9')
     

Notice how the hash code did not change because the force field definition has not physically changed. In contrast, if we change a physical property such as the `defs`, its hash code changes i.e.

.. code-block:: python

    >>> ff = mmelemental.models.ForceField(name="my_forcefield", defs=["OW"])
    >>> ff
     ForceField(name='my_forcefield', form=[], hash='0d9fc00')

Here we define an entity called "OW" which refers to an atom type of Oxygen.

Functional form
===============
Here we discuss the use of the ForceField_ class for constructing force fields such as Amber99 and Sage OFF.
To be completed.


Assigned parameters
===================
ForceField_ objects can also store force field parameters that were assigned to a particular molecule.
To be completed.


I/O operations
--------------
File operations
===============
Coming soon.

Data conversion
===============
Coming soon.

Coarse-graining
---------------
The ForceField_ models are applicable to any kind of particle systems i.e. the underlying object these models describe does not have to be atomic. The `symbols` property could
for instance represent entities rather than atoms (although this will negate atomic properties such as atomic or mass numbers). 

Validation
----------
MMElemental performs only data type validation on any constructed model. However, beyond basic validation and sanity checks, MMElemental does not perform any scientific validation. This is what enables 
MMElemental to support coarse-graining for instance. For domain-specific (i.e. scientific) validation, MMElemental can theoretically make use of MMIC_ validators similarly to how it uses translators to 
parse and write to various file formats.
