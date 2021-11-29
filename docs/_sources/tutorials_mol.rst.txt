Structure
#########
.. _MMSchema: https://molssi.github.io/mmschema
.. _pydantic: https://sphinx-pydantic.readthedocs.io
.. _Molecule: molecules.html#molecule
.. _Topology: molecules.html#topology
.. _Jupyter: https://jupyter.org
.. _nglview: http://nglviewer.org/nglview/latest/
.. _pint: https://pint.readthedocs.io
.. _MMIC: https://github.com/MolSSI/mmic

The most fundamental and general model for representing molecules in MMElemental is with the Molecule_ model. 

Overview & Instantiation
------------------------

Molecule_ objects can be created by supplying keyword arguments to the constructor. Per the MMSchema_ specification, the only required keyword argument is `symbols` while all other fields
being optional. Hence, the simplest way to create a molecule is via:

.. code-block:: python

    >>> mol = mmelemental.models.Molecule(symbols = ["O", "H", "H"])

This creates a Molecule_ object with a unique hash code and a name that defaults to the molecular formula (can be overwritten with the keyword arg `name`).

.. code-block:: python

    >>> mol
     Molecule(name='H2O', hash='83a5d8a')

We can visualize this molecule if we're using Jupyter_ notebook/lab with nglview_. In order to do that, we need to specify the atomic positions by passing the `geometry` keyword to the constructor i.e.

.. code-block:: python

    >>> mol = mmelemental.models.Molecule(
    >>>     symbols = ["O", "H", "H"],
    >>>     geometry = [2.0, 2.09, 0.0, 2.82, 2.09, 0.58, 1.18, 2.09, 0.58],
    >>> )

    >>> mol
     Molecule(name='H2O', hash='631741f')

Notice how the hash code changed because we've specified another core field (`geometry`) in the molecule model i.e. the molecule definition has physically changed. In contrast, if we change a property such as
the `name` of the molecule, its hash code remains the same i.e.

.. code-block:: python

    >>> mol = mmelemental.models.Molecule(
    >>>     name = "water",
    >>>     symbols = ["O", "H", "H"],
    >>>     geometry = [2.0, 2.09, 0.0, 2.82, 2.09, 0.58, 1.18, 2.09, 0.58],
    >>> )

    >>> mol
     Molecule(name='H2O', hash='631741f')

We can visualize this molecule if nglview_ is installed:

.. code-block:: python

    >>> mol.show()

or

    >>> import nglview
    >>> nglview.show_mmelemental(mol)

NGLview guesses the connectivity or bonds in a given molecule. We can, however, specify the connectivity (and the bond order).

.. code-block:: python

    >>> mol = mmelemental.models.Molecule(
    >>>     symbols = ["O", "H", "H"],
    >>>     geometry = [2.0, 2.09, 0.0, 2.82, 2.09, 0.58, 1.18, 2.09, 0.58],
    >>>     connectivity = [(0, 1, 1.0), (0, 2, 1.0)],
    >>> )

    >>> mol
     Molecule(name='H2O', hash='18705f4')

Since connectivity changes the physical definition of a molecule, the hash code changes when connectivity is specified or modified. Note that the `connectivity` field should always
be a list of tuple of the form (atom1_index, atom2_index, bond_order). Otherwise, pydantic_ will throw in a validation error.

Under the hood, all array fields are cast as numpy arrays, and every field set by the user becomes accessible (but cannot be modified since `Molecule` is an immutable object). For instance,
we can access the connectivity or geometry via:

.. code-block:: python

    >>> mol.geometry
     array([2.  , 2.09, 0.  , 2.82, 2.09, 0.58, 1.18, 2.09, 0.58])

    >>> mol.connectivity
     array([(0, 1, 1.), (0, 2, 1.)], dtype=[('f0', '<i8'), ('f1', '<i8'), ('f2', '<f8')])    

We can also access other attributes we did not explicitly specify such as atomic numbers and dimensionality:

.. code-block:: python

    >>> mol.atomic_numbers
     array([8, 1, 1])

    >>> mol.ndim
     3

Certain physical properties such as `geometry`, `velocities`, `masses`, and `molecular_charge` have default units fields as well that can be set based on physically consistent and supported units
available in pint_. For example, we can access the defaul geometry unit (`angstrom`) or access all default units available in this model:

.. code-block:: python

    >>> mol.geometry_units
     'angstrom'     

    >>> mmelemental.models.Molecule.default_units
     {
       'masses_units': 'unified_atomic_mass_unit', 'molecular_charge_units': 'elementary_charge', 
       'formal_charges_units': 'elementary_charge', 'partial_charges_units': 'elementary_charge', 
       'geometry_units': 'angstrom', 'velocities_units': 'angstrom / femtosecond'
     }

The instance method `mol.units` returns all stored units that were assigned to the object. For example:

.. code-block:: python

    >>> mol = mmelemental.models.Molecule(
    >>>     symbols = ["O", "H", "H"],
    >>>     geometry = [0.2, 0.209, 0.0, 0.282, 0.209, 0.58, 0.118, 0.209, 0.58],
    >>>     geometry_units = "nm",
    >>>     connectivity = [(0, 1, 1.0), (0, 2, 1.0)],
    >>> )

    >>> mol
     Molecule(name='H2O', hash='720bf3a')

In this case, we specified the geometry in nanometers, which is not the default geometry unit. Note how the hash code changed as well despite the fact that the molecule remains scientifically 
(but not programmatically) equivalent. Also note that `mol.units != mol.default_units` since the former has its geometry units in nanometers whereas `mol.default_units` is a function of only
the model schema itself (and not any instant of the class).

In the next section, we will go over how molecules can be created from files, seralized, and written to files.

Topological data
----------------

A Topology_ model in MMElemental is a subset of Molecule_ which captures abstract information about a molecule's connectivity. A Topology_ object
does not, however, store particle positions. This model is particularly useful for efficiently storing trajectories generated from classical molecular dynamics, in which
the connectivity of a molecule is constant while the particle positions change over time. We can create a Topology_ object from an existing molecule with the `get_topology` 
method as shown:

.. code-block:: python

    >>> top = mol.get_topology()

    >>> top
     Topology(name='top_from_mol', hash='6915dc7')

Alternatively, 

A Topology_ object can also be used to instantiate a Molecule_ object. For instance,

.. code-block:: python

    >>> mmelemental.models.Molecule(**top.dict(exclude={"schema_name"}))
     Molecule(name='top_from_mol', hash='a5f83e3')

Notice that `top.dict(exclude={"schema_name"})` extracts all populated fields and returns them in a python dictionary, excluding the `schema_name` (which is by default `mmschema_molecule` for molecules).


I/O operations
--------------
File operations
===============
To be completed.

Data conversion
===============
To be completed.

Coarse-graining
---------------
The Molecule_ and Topology_ models are applicable to any kind of particle systems i.e. the underlying object these models describe does not have to be atomic. The `symbols` property could
for instance represent entities rather than atoms (although this will negate atomic properties such as atomic or mass numbers). 

Validation
----------
MMElemental performs data type validation on any constructed model. However, beyond basic validation and sanity checks, MMElemental does not perform any scientific validation. This is why non-atomic entities
are supported for instance for coarse-graining, and this is what enables MMElemental to support a wide range of applications. For domain-specific (i.e. scientific) validation, MMElemental can theoretically 
make use of MMIC_ validators similarly to how it uses translators to parse/write various file formats.
