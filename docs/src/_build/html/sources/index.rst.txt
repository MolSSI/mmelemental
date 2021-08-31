.. mmelemental documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _MMSchema: https://molssi.github.io/mmschema
.. _pydantic: https://sphinx-pydantic.readthedocs.io

MMElemental documentation
=========================

MMElemental is an implementation of MMSchema_. Specifically, MMElemental is a python package that provides pydantic_ models for molecular mechanics (MM).
MMElemental makes it easy to standardize MM applications and workflows in python based on the MMSchema specification, and it provides interoperability capabilities 
that allow conversion between different representations of MM objects or file formats.


.. image:: _static/mmel_mmschema.png
   :scale: 55 %
   :align: center

.. toctree::
   :maxdepth: 2
   :caption: Introduction:

   get_started
   models
   components

.. toctree::
   :maxdepth: 1
   :caption: Modules:

   molecules
   forcefields
   collections

.. toctree::
   :maxdepth: 1
   :caption: Tutorials:

   tutorials_mol
   tutorials_ff
   tutorials_traj
   tutorials_debug

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
