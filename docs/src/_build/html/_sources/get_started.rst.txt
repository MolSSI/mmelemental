Prerequisites
#############

.. _pydantic: https://sphinx-pydantic.readthedocs.io
.. _PyPi: https://pypi.org

Quick Installation
==================
MMElemental will be available from PyPi_ i.e. you can easily install it with `pip`.

.. code-block:: bash

   pip install mmelemental

Dev Installation
================
For developers who would like to install MMElemental from source::

   pip install git+https://github.com/MolSSI/mmelemental

Requirements
============
The core and optional basic requirements in MMElemental are:

- pydantic_ (required)
- `cmselemental <https://pypi.org/project/cmselemental/>`_ (required)
- `qcelemental <https://pypi.org/project/qcelemental/>`_ (optional)
- `pytest <https://pytest.org>`_ (optional)
- `mm_data <https://github.com/MolSSI/mm_data>`_ (optional)

Because of the modular design in MMElemental, the number of optional requirements can be infinite. For instance, parsing a specific file format such as mmCIF or PDB
requires having mmic_translator installed (and a tactic component thereof that can parse mmCIF or PDB files). MMElemental informs the user in runtime whenever 
a particular package is needed for a requested task.
