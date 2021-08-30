Molecule
########
.. _MMSchema: https://molssi.github.io/mmschema
.. _pydantic: https://sphinx-pydantic.readthedocs.io
.. _Molecule: molecules.html#molecule
.. _Topology: molecules.html#topology
.. _Jupyter: https://jupyter.org
.. _nglview: http://nglviewer.org/nglview/latest/
.. _pint: https://pint.readthedocs.io


The most fundamental and general model for representing molecules in MMElemental is with the Molecule_ model. 

Overview & Instantiation
------------------------

Molecule_ objects can be created by supplying keyword arguments to the constructor. Per the MMSchema_ specification, the only required keyword argument is `symbols` with all other fields
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
available in pint_. For example, we can access the defaul geometry unit (`angstrom`) or access all units available in this model.

.. code-block:: python

    >>> mol.geometry_units
     'angstrom'     

    >>> mol.get_units()
     {'masses_units': 'amu', 'molecular_charge_units': 'e', 'geometry_units': 'angstrom', 'velocities_units': 'angstrom/fs'}

In the next section, we will go over how molecules can be created from files, seralized and written to files.

Topological data
----------------

A Topology_ model in MMElemental is a subset of Molecule_ which captures symbolic information such as particle symbols and connectivity.
We can create a Topology_ object from an existing molecule with the `get_topology` method as shown:

.. code-block:: python

    >>> top = mol.get_topology()

    >>> top
     Topology(name='top_from_mol', hash='6915dc7')

Alternatively, 

A Topology_ object can also be used to instantiate a Molecule_ object. For instance,

.. code-block:: python

    >>> mmelemental.models.Molecule(**top.dict(exclude={"schema_name"}))
     Molecule(name='top_from_mol', hash='a5f83e3')

Notice that `top.dict(exclude={"schema_name"})` extracts all populated fields and return them in a python dictionary, excluding the `schema_name` (which is by default `mmschema_molecule` for molecules).


I/O operations
--------------
File operations
===============

Data conversion
===============

Coarse-graining
---------------

