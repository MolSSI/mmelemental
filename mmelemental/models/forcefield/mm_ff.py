from pydantic import Field
import importlib
from typing import Any, List, Optional
from qcelemental.models.types import Array

# Import MM models
from mmelemental.models.base import ProtoModel
from mmelemental.models.base import ToolkitModel
from mmelemental.models.util.output import FileOutput

# Import MM components
from mmelemental.components.trans import TransComponent

__all__ = ["ForceField"]


class Bonds(ProtoModel):
    length: Optional[Array[float]] = Field(
        None, description="Bond equilibrium lengths. Default unit is Angstroms."
    )
    length_units: Optional[str] = Field("angstroms", description="Bond lengths unit.")
    spring: Optional[Array[float]] = Field(
        0, description="Bond spring constant. Default unit is kJ/(mol*angstrom**2)."
    )
    spring_units: Optional[str] = Field(
        "kJ/(mol*angstrom**2)", description="Bond spring constant unit."
    )
    params: Optional[Array[float]] = Field(
        None,
        description="Extra or custom parameters for describing the bond potential.",
    )
    form: Optional[str] = Field(
        None, description="Bond potential form e.g. harmonic, morse, etc."
    )


class Angles(ProtoModel):
    angle: Optional[Array[float]] = Field(
        None, description="Equilibrium angles. Default unit is degrees."
    )
    angle_units: Optional[str] = Field(
        "degrees", description="Equilibrium angle units."
    )
    spring: Optional[Array[float]] = Field(0, description="Angle spring constant. ")
    spring_units: Optional[str] = Field(
        "kJ/(mol*degrees**2)", description="Angle spring constant unit."
    )
    params: Optional[Array[float]] = Field(
        None,
        description="Extra or custom parameters for describing the angle potential.",
    )
    form: Optional[str] = Field(
        None, description="Angle potential form e.g. harmonic, quartic, etc."
    )


class Dihedrals(ProtoModel):
    dihedrals: Optional[Array[Array[float]]] = Field(
        None,
        description="Dihedral/torsion parameters: eq energy, phase in degrees, extra params.",
    )
    dihedrals_type: Optional[List[str]] = Field(
        None,
        description="Proper dihedral potential form e.g. harmonic, helix, fourier, etc.",
    )


class ImproperDihedrals(ProtoModel):
    improper_dihedrals: Optional[Array[Array[float]]] = Field(
        None, description="Improper dihedral/torsion parameters."
    )
    improper_dihedrals_type: Optional[List[str]] = Field(
        None,
        description="Improper dihedral potential form e.g. harmonic, fourier, etc.",
    )


class NonBonded(ProtoModel):
    epsilon: Optional[Array[float]] = Field(
        0,
        description="The epsilon (well depth) Lennard-Jones parameter. Default unit is kJ/mol.",
    )
    epsilon_units: Optional[str] = Field(
        "kJ/mol", description="The epsilon (well depth) Lennard-Jones unit."
    )
    sigma: Optional[Array[float]] = Field(
        0,
        description="The distance at which the Lennard-Jones potential is 0. Default unit is angstroms.",
    )
    sigma_units: Optional[str] = Field(
        "angstrom", description="The Lennard-Jones sigma unit."
    )
    params: Optional[Array[float]] = Field(
        None, description="Extra or custom non-bonded short potential parameters."
    )
    form: Optional[str] = Field(
        None,
        description="Non-bonded short potential form e.g. Lennard-Jones, Buckingham, etc.",
    )
    charges: Optional[Array[float]] = Field(
        None, description="Atomic charges. Default unit is in elementary charge units."
    )
    charges_units: Optional[str] = Field("e", description="Atomic charge unit.")


class ForceField(ProtoModel):
    bonds: Optional[Bonds] = Field(None, description="2-body covalent bond schema.")
    angles: Optional[Angles] = Field(None, description="3-body angular bond schema.")
    dihedrals: Optional[Dihedrals] = Field(
        None, description="4-body torsional bond schema."
    )
    im_dihedrals: Optional[ImproperDihedrals] = Field(
        None, description="Improper dihedral bond schema."
    )
    nonbonded: NonBonded = Field(..., description="Non-bonded parameters schema.")
    charge_groups: Optional[Array[int]] = Field(
        None, description="Charge groups per atom. Length of the array must be natoms."
    )
    exclusions: Optional[str] = Field(
        None,
        description="Which pairs of bonded atoms to exclude from non-bonded calculations. \
    	The rules to apply in choosing bonded exclusions are specifed in the configuration file using the exclude parameter. The \
    	choices for exclusions are None, 1-2, 1-3, 1-4, etc. With None, no atom pairs are excluded. With 1-2, only atoms that are connected \
    	via a linear bond are excluded. With 1-3, any pair of atoms connected via a bond or bonded to a common third atom are excluded. \
    	With 1-4, any atoms that are connected by a linear bond, or a sequence of two bonds, or a sequence of three bonds, are excluded. \
    	With scaled1-4, exclusions are applied in the same manner as the 1-3 setting, but those pairs that are connected by a sequence of \
    	3 bonds are calculated using the modified 1-4 methods described rather than the standard force calculations.",
    )
    inclusions: Optional[str] = Field(
        None,
        description="Which pairs of 1-4 excluded bonded atoms to include in non-bonded calculations.",
    )
    name: Optional[str] = Field(
        None, description="Forcefield name e.g. charmm27, amber99, etc."
    )
    types: Optional[List[str]] = Field(
        None,
        description="Atom types e.g. HH31. The type names are associated with the atomic \
        elements defined in other objects e.g. see the :class:``Molecule`` model.",
    )

    @classmethod
    def from_file(
        cls, filename: str, dtype: Optional[str] = None, **kwargs
    ) -> "ForceField":
        """
        Constructs a ForceField object from a file.
        Parameters
        ----------
        filename : str
            The topology or FF filename to build from.
        dtype : str, optional
            The type of file to interpret. If not set, mmelemental attempts to discover the file type.
        **kwargs
            Any additional keywords to pass to the constructor.
        Returns
        -------
        ForceField
            A constructed ForceField object.
        """
        tkmol = cls._tkFromFile(filename=filename, dtype=dtype)
        return cls.from_data(tkmol, dtype=tkmol.dtype)

    @classmethod
    def from_data(cls, data: Any, **kwargs) -> "ForceField":
        """
        Constructs a ForceField object from a data object.
        Parameters
        ----------
        data: Any
            Data to construct ForceField from.
        **kwargs
            Additional kwargs to pass to the constructors.
        Returns
        -------
        ForceField
            A constructed ForceField object.
        """
        return data.to_schema(**kwargs)

    def to_file(self, filename: str, dtype: Optional[str] = None, **kwargs) -> None:
        """Writes the ForceField to a file.
        Parameters
        ----------
        filename : str
            The filename to write to
        dtype : Optional[str], optional
            The type of file to write (e.g. psf, top, etc.), attempts to infer dtype from
            file extension if not provided.
        **kwargs
            Additional kwargs to pass to the constructor.
        """
        if not dtype:
            from pathlib import Path

            ext = Path(filename).suffix

        translator = TransComponent.find_ffwrite_tk(ext)
        tkff = self._tkFromSchema(translator=translator, **kwargs)
        tkff.to_file(filename, dtype=ext.strip("."), **kwargs)  # pass dtype?

    def to_data(self, dtype: str, **kwargs) -> ToolkitModel:
        """Converts ForceField to toolkit-specific forcefield object.
        Parameters
        ----------
        dtype: str
            The type of data object to convert to e.g. parmed, MDAnalysis, etc.
        **kwargs
            Additional kwargs to pass to the constructor.
        Returns
        -------
        ToolkitModel
            Toolkit-specific forcefield model
        """
        return self._tkFromSchema(dtype=dtype, **kwargs)

    def _tkFromSchema(
        self, translator: str = None, dtype: str = None, **kwargs
    ) -> ToolkitModel:
        """
        Helper function that constructs a toolkit-specific forcefield from MMSchema ForceField.
        Which toolkit-specific component is called depends on which package is installed on the system.
        Parameters
        ----------
        translator: str, optional
            Translator name e.g. mmic_mda, etc.
        dtype: str, optional
            Data type e.g. MDAnalysis, parmed, etc.
        **kwargs
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

    @classmethod
    def _tkFromFile(cls, filename: str, dtype: str = None, **kwargs) -> ToolkitModel:
        """Helper function that constructs a toolkit-specific parameterized molecule from an input file.
        Which toolkit-specific component is called depends on which package is installed on the system.
        Parameters
        ----------
        filename: str
            Topology file name. The file should define FF parameters.
        dtype: str, optional
            The type of data file object to read e.g. top, psf, etc.
        **kwargs
            Additional kwargs to pass to the constructor.
        Returns
        -------
        ToolkitModel
        """

        fileobj = FileOutput(path=filename)
        dtype = dtype or fileobj.ext
        translator = TransComponent.find_ffread_tk(dtype)

        if not translator:
            raise ValueError(
                f"Could not read file with ext {dtype}. Please install an appropriate TransComponent."
            )

        if importlib.util.find_spec(translator):
            mod = importlib.import_module(translator)
            tkff = mod._classes_map.get("ForceField")

        if not tkff:
            raise ValueError(
                f"No ForceField model found while looking in translator: {translator}."
            )

        return tkff.from_file(filename=fileobj.abs_path, dtype=dtype.strip("."))
