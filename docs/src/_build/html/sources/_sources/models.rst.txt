Models
######
.. _MMSchema: https://molssi.github.io/mmschema
.. _pydantic: https://sphinx-pydantic.readthedocs.io

MMElemental provides models based on the MMSchema_ specification. These models are immutable,
and they provide serialization and data validation based on the pydantic_ library. Each model
has a set of fields that, when suitable, are used to automatically generate a unique hash 
that enables each object to be uniquely identified.
