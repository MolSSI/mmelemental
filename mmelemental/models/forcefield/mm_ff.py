from pydantic import Field
from typing import Any, Tuple, Union, List, Dict, Optional
from mmelemental.models.base import Base
from qcelemental.models.types import Array
from mmelemental.models.util.input import FileInput
from mmic.components.blueprints.generic_component import GenericComponent
from .io_ff import FFInput


class ForceField(Base):
    bonds: Optional[Array[Array[float]]] = Field(
        None,
        description="Bond parameters: eq length in Angstrom, constant (e.g. spring const) in kJ/mol*Angstrom^2, extra params.",
    )
    bonds_type: Optional[List[str]] = Field(
        None, description="Bond potential form e.g. harmonic, morse, etc."
    )

    angles: Optional[Array[Array[float]]] = Field(
        None,
        description="Angle parameters: eq angle in degrees, constant (e.g. spring const) in kJ/mol*Angstrom, extra params.",
    )
    angles_type: Optional[List[str]] = Field(
        None, description="Angle potential form e.g. harmonic, quartic, etc."
    )

    dihedrals: Optional[Array[Array[float]]] = Field(
        None,
        description="Dihedral/torsion parameters: eq energy in kJ/mol, phase in degrees, extra params.",
    )
    dihedrals_type: Optional[List[str]] = Field(
        None,
        description="Proper dihedral potential form e.g. harmonic, helix, fourier, etc.",
    )

    improper_dihedrals: Optional[Array[Array[float]]] = Field(
        ..., description="Dihedral/torsion parameters."
    )
    improper_dihedrals_type: Optional[List[str]] = Field(
        None,
        description="Improper dihedral potential form e.g. harmonic, fourier, etc.",
    )

    charges: Array[float] = Field(
        ..., description="Atomic charges in elementary charge units."
    )
    charge_groups: Optional[Array[int]] = Field(
        None, description="Charge groups per atom."
    )

    nonbonded: Array[float] = Field(
        ..., description="Non-bonded short potential parameters."
    )
    nonbonded_type: Optional[List[str]] = Field(
        None,
        description="Non-bonded short potential form e.g. Lennard-Jones, Buckingham, etc.",
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
        cls, filename: Union[FileInput, str], dtype: Optional[str] = None, **kwargs
    ) -> "ForceField":
        """
        Constructs a ForceField object from a file.
        Parameters
        ----------
        filename : str
            The filename to build from
        dtype : str, optional
            The type of file to interpret. If not set, mmelemental attempts to discover the file type.
        **kwargs
            Any additional keywords to pass to the constructor
        Returns
        -------
        ForceField
            A constructed ForceField class.
        """
        if not isinstance(filename, FileInput):
            filename = FileInput(path=filename)

        ff_input = FFInput(file=filename)
        return FFReaderComponent.compute(ff_input)

    @classmethod
    def from_data(cls, data: Any, **kwargs: Dict[str, Any]) -> "ForceField":
        """
        Constructs a ForceField object from a data object.
        Parameters
        ----------
        data: Any
            Data to construct ForceField from
        **kwargs
            Additional kwargs to pass to the constructors. kwargs take precedence over data.
        Returns
        -------
        ForceField
            A constructed ForceField class.
        """
        ff_input = FFInput(data=data)
        return FFReaderComponent.compute(ff_input)


class FFReaderComponent(GenericComponent):
    """Factory component that constructs a Molecule object from MolInput.
    Which toolkit-specific component is called depends on data type and
    which toolkits are installed on the system."""

    @classmethod
    def input(cls):
        return FFInput

    @classmethod
    def output(cls):
        return ForceField

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, Dict[str, Any]]:

        if inputs.data:
            dtype = inputs.data.dtype
            if dtype == "parmed":
                from mmelemental.components.trans.parmed_component import (
                    ParmedToFFComponent,
                )

                return True, ParmedToFFComponent.compute(inputs)
            else:
                raise NotImplementedError(f"Data type not yet supported: {dtype}.")
        elif inputs.file:
            toolkit = inputs.files_toolkit()
            if toolkit == "parmed":
                from mmelemental.components.trans.parmed_component import (
                    ParmedToFFComponent,
                )

                return True, ParmedToFFComponent.compute(inputs)
        else:
            raise NotImplementedError(
                "Forcefield can be instantiated only from files or other ssupported data objects."
            )
