from pydantic import Field, constr
import importlib
import hashlib
import json
from typing import Any, List, Dict, Optional
import qcelemental
import numpy

# MM models
from mmelemental.models.base import ProtoModel, Provenance, provenance_stamp
from mmelemental.models.util.output import FileOutput
from .nonbonded import NonBonded
from .bonded import Bonds, Angles

# Generic translator component
from mmic_translator.models.base import ToolkitModel
from mmic_translator.components import TransComponent


mmschema_forcefield_default = "mmschema_forcefield"


__all__ = ["ForceField"]


class Dihedrals(ProtoModel):
    angle: Optional[qcelemental.models.types.Array[float]] = Field(
        None, description="Equilibrium angles. Default unit is degrees."
    )
    angle_units: Optional[str] = Field(
        "degrees", description="Equilibrium angle units."
    )
    spring: Optional[qcelemental.models.types.Array[float]] = Field(
        0, description="Dihedral spring constant. "
    )
    spring_units: Optional[str] = Field(
        "kJ/(mol*degrees**2)", description="Dihedral spring constant unit."
    )
    params: Optional[qcelemental.models.types.Array[float]] = Field(
        None,
        description="Extra or custom parameters for describing the dihedral potential.",
    )
    form: Optional[str] = Field(
        None, description="Dihedral potential form e.g. harmonic, fourier, etc."
    )


class ImproperDihedrals(ProtoModel):
    im_dihedrals: Optional[
        qcelemental.models.types.Array[qcelemental.models.types.Array[float]]
    ] = Field(None, description="Improper dihedral/torsion parameters.")
    im_dihedrals_type: Optional[List[str]] = Field(
        None,
        description="Improper dihedral potential form e.g. harmonic, fourier, etc.",
    )


class ForceField(ProtoModel):
    schema_name: constr(
        strip_whitespace=True, regex="^(mmschema_forcefield)$"
    ) = Field(  # type: ignore
        mmschema_forcefield_default,
        description=(
            f"The MMSchema specification to which this model conforms. Explicitly fixed as {mmschema_forcefield_default}."
        ),
    )
    schema_version: int = Field(  # type: ignore
        0,
        description="The version number of ``schema_name`` to which this model conforms.",
    )
    bonds: Optional[Bonds] = Field(  # type: ignore
        None, description="2-body covalent bond model."
    )
    angles: Optional[Angles] = Field(  # type: ignore
        None, description="3-body angular bond model."
    )
    dihedrals: Optional[Dihedrals] = Field(  # type: ignore
        None, description="4-body torsional bond model."
    )
    im_dihedrals: Optional[ImproperDihedrals] = Field(  # type: ignore
        None, description="Improper dihedral bond model."
    )
    nonbonded: NonBonded = Field(  # type: ignore
        ..., description="Non-bonded parameters model."
    )
    charge_groups: Optional[qcelemental.models.types.Array[int]] = Field(
        None, description="Charge groups per atom. Length of the array must be natoms."
    )
    exclusions: Optional[str] = Field(  # type: ignore
        None,
        description="Which pairs of bonded atoms to exclude from non-bonded calculations. \
    	The rules to apply in choosing bonded exclusions are specifed in the configuration file using the exclude parameter. The \
    	choices for exclusions are None, 1-2, 1-3, 1-4, etc. With None, no atom pairs are excluded. With 1-2, only atoms that are connected \
    	via a linear bond are excluded. With 1-3, any pair of atoms connected via a bond or bonded to a common third atom are excluded. \
    	With 1-4, any atoms that are connected by a linear bond, or a sequence of two bonds, or a sequence of three bonds, are excluded. \
    	With scaled1-4, exclusions are applied in the same manner as the 1-3 setting, but those pairs that are connected by a sequence of \
    	3 bonds are calculated using the modified 1-4 methods described rather than the standard force calculations.",
    )
    inclusions: Optional[str] = Field(  # type: ignore
        None,
        description="Which pairs of 1-4 excluded bonded atoms to include in non-bonded calculations.",
    )
    name: Optional[str] = Field(  # type: ignore
        None, description="Forcefield name e.g. charmm27, amber99, etc."
    )
    types: Optional[List[str]] = Field(  # type: ignore
        None,
        description="Atom types e.g. HH31. The type names are associated with the atomic \
        elements defined in other objects e.g. see the :class:``Molecule`` model.",
    )
    symbols: Optional[List[str]] = Field(  # type: ignore
        None, description="An ordered (natom,) list of atomic elemental symbols."
    )
    masses_: Optional[List[float]] = Field(  # type: ignore
        None,
        description="List of atomic masses. If not provided, the mass of each atom is inferred from its most common isotope. "
        "If this is provided, it must be the same length as ``symbols``.",
    )
    masses_units: Optional[str] = Field(  # type: ignore
        "amu",
        description="Units for atomic masses. Defaults to unified atomic mass unit.",
    )
    atomic_numbers_: Optional[
        qcelemental.models.types.Array[numpy.int16]
    ] = Field(  # type: ignore
        None,
        description="An optional ordered 1-D array-like object of atomic numbers of shape (nat,). Index "
        "matches the 0-indexed indices of all other per-atom settings like ``symbols`` and ``real``. "
        "Values are inferred from the ``symbols`` list if not explicitly set. "
        "Ghostedness should be indicated through ``real`` field, not zeros here.",
        shape=["nat"],
    )
    # Extras
    provenance: Provenance = Field(  # type: ignore
        provenance_stamp(__name__),
        description="The provenance information about how this object (and its attributes) were generated, "
        "provided, and manipulated.",
    )
    extras: Dict[str, Any] = Field(  # type: ignore
        None,
        description="Additional information to bundle with the molecule. Use for schema development and scratch space.",
    )
    class Config(ProtoModel.Config):
        serialize_skip_defaults = True
        repr_style = lambda self: [("name", self.name), ("hash", self.get_hash()[:7])]
        fields = {"masses_": "masses", "atomic_numbers_": "atomic_numbers"}

        def schema_extra(schema, model):
            # below addresses the draft-04 issue until https://github.com/samuelcolvin/pydantic/issues/1478 .
            schema["$schema"] = "http://json-schema.org/draft-04/schema#"

    # Properties
    @property
    def masses(self) -> qcelemental.models.types.Array[float]:
        masses = self.__dict__.get("masses_")
        if masses is None:
            masses = numpy.array(
                [qcelemental.periodictable.to_mass(x) for x in self.symbols]
            )
        return masses

    @property
    def atomic_numbers(self) -> qcelemental.models.types.Array[numpy.int16]:
        atomic_numbers = self.__dict__.get("atomic_numbers_")
        if atomic_numbers is None:
            atomic_numbers = numpy.array(
                [qcelemental.periodictable.to_Z(x) for x in self.symbols]
            )
        return atomic_numbers

    # Constructors
    @classmethod
    def from_file(
        cls,
        filename: str,
        dtype: Optional[str] = None,
        translator: Optional[str] = None,
        **kwargs,
    ) -> "ForceField":
        """
        Constructs a ForceField object from a file.
        Parameters
        ----------
        filename: str
            The topology or FF filename to build from.
        dtype: Optional[str], optional
            The type of file to interpret e.g. psf. If unset, mmelemental attempts to discover the file type.
        translator: Optional[str], optional
            Translator name e.g. mmic_parmed. Takes precedence over dtype. If unset, MMElemental attempts 
            to find an appropriate translator if it is registered in the :class:``TransComponent`` class. 
        **kwargs: Optional[Dict[str, Any]], optional
            Any additional keywords to pass to the constructor.
        Returns
        -------
        ForceField
            A constructed ForceField object.
        """

        fileobj = FileOutput(path=filename)
        dtype = dtype or fileobj.ext.strip(".")
        ext = "." + dtype

        if not translator:
            translator = TransComponent.find_ffread_tk(ext)

        if not translator:
            raise ValueError(
                f"Could not read file with ext {dtype}. Please install an appropriate TransComponent."
            )

        if importlib.util.find_spec(translator):
            mod = importlib.import_module(translator)
            tkff_class = mod._classes_map.get("ForceField")

        if not tkff_class:
            raise ValueError(
                f"No ForceField model found while looking in translator: {translator}."
            )

        tkff = tkff_class.from_file(filename=fileobj.abs_path, dtype=dtype)

        return cls.from_data(tkff, dtype=tkff.dtype)

    @classmethod
    def from_data(cls, data: Any, **kwargs) -> "ForceField":
        """
        Constructs a ForceField object from a data object.
        Parameters
        ----------
        data: Any
            Data to construct ForceField from.
        **kwargs: Optional[Dict[str, Any]], optional
            Additional kwargs to pass to the constructors.
        Returns
        -------
        ForceField
            A constructed ForceField object.
        """
        return data.to_schema(**kwargs)

    def to_file(
        self,
        filename: str,
        dtype: Optional[str] = None,
        translator: Optional[str] = None,
        **kwargs: Dict[str, Any],
    ) -> None:
        """Writes the ForceField to a file.
        Parameters
        ----------
        filename : str
            The filename to write to
        dtype : Optional[str], optional
            The type of file to write (e.g. psf, top, etc.), attempts to infer dtype from
            file extension if not provided.
        translator: Optional[str], optional
            Translator name e.g. mmic_parmed. Takes precedence over dtype. If unset, MMElemental attempts 
            to find an appropriate translator if it is registered in the :class:``TransComponent`` class. 
        **kwargs: Optional[str, Dict], optional
            Additional kwargs to pass to the constructor.
        """
        if not dtype:
            from pathlib import Path

            ext = Path(filename).suffix
        else:
            ext = "." + dtype

        mode = kwargs.pop("mode", "w")

        if ext == ".json":
            stringified = self.json(**kwargs)

            with open(filename, mode) as fp:
                fp.write(stringified)

            return

        if not translator:
            translator = TransComponent.find_ffwrite_tk(ext)

        if not translator:
            raise NotImplementedError(
                f"File extension {ext} not supported with any installed translators."
            )

        tkff = self.to_data(translator=translator, **kwargs)
        tkff.to_file(filename, dtype=dtype, **kwargs)  # pass dtype?

    def to_data(
        self,
        dtype: Optional[str] = None,
        translator: Optional[str] = None,
        **kwargs: Dict[str, Any],
    ) -> ToolkitModel:
        """
        Constructs a toolkit-specific forcefield from MMSchema ForceField.
        Which toolkit-specific component is called depends on which package is installed on the system.
        Parameters
        ----------
        translator: Optional[str], optional
            Translator name e.g. mmic_parmed. Takes precedence over dtype. If unset, MMElemental attempts 
            to find an appropriate translator if it is registered in the :class:``TransComponent`` class. 
        dtype: Optional[str], optional
            Data type e.g. MDAnalysis, parmed, etc.
        **kwargs: Optional[Dict[str, Any]]
            Additional kwargs to pass to the constructors.
        Results
        -------
        ToolkitModel
            Toolkit-specific ForceField object
        """

        if not translator:
            if not dtype:
                raise ValueError(
                    f"Either translator or dtype must be supplied when calling {__name__}."
                )
            translator = TransComponent.find_trans(dtype)

        if importlib.util.find_spec(translator):
            mod = importlib.import_module(translator)
            tkff = mod._classes_map.get("ForceField")

            if not tkff:
                raise ValueError(
                    f"No ForceField model found while looking in translator: {translator}."
                )

            return tkff.from_schema(self)
        else:
            raise NotImplementedError(
                f"Translator {translator} not available. Make sure it is properly installed."
            )

    def __eq__(self, other):
        """
        Checks if two molecules are identical. This is a molecular identity defined
        by scientific terms, and not programing terms, so it's less rigorous than
        a programmatic equality or a memory equivalent `is`.
        """

        if isinstance(other, dict):
            other = Forcefield(**other)
        elif isinstance(other, Forcefield):
            pass
        else:
            raise TypeError(
                f"Comparison molecule not understood of type '{type(other)}'."
            )

        return self.get_hash() == other.get_hash()

    @property
    def hash_fields(self):
        return [
            "nonbonded",
            "bonds",
            "angles",
            "dihedrals",
            "im_dihedrals",
            "exclusions",
            "inclusions",
        ]

    def get_hash(self):
        """
        Returns the hash of the force field object.
        """

        m = hashlib.sha1()
        concat = ""

        # np.set_printoptions(precision=16)
        for field in self.hash_fields:
            data = getattr(self, field)
            if data is not None:
                # if field == "nonbonded":
                #    data = qcelemental.models.molecule.float_prep(data, GEOMETRY_NOISE)

                concat += json.dumps(data, default=lambda x: x.ravel().tolist())

        m.update(concat.encode("utf-8"))
        return m.hexdigest()
