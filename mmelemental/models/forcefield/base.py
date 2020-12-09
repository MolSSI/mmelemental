from pydantic import Field
from typing import Any, Tuple, Union, List, Dict, Optional
from mmelemental.models.base import Base
from qcelemental.models.types import Array

class ForceField(Base):
    bonds: Optional[Tuple[int, int, float]] = Field(None, description = 'Bond parameters.')
    bonds_pot: Optional[str] = Field(None, description = 'Bond potential form e.g. harmonic, morse, etc.')

    angles: Optional[Tuple[int, int, float]] = Field(None, description = 'Angle parameters.')
    angles_pot: Optional[str] = Field(None, description = 'Angle potential form e.g. harmonic, quartic, etc.')

    proper_dihedrals: Optional[Tuple[int, int, int, int, float]] = Field(None, description = 'Dihedral/torsion parameters.')
    proper_dihedrals_pot: Optional[str] = Field(None, description = 'Proper dihedral potential form e.g. harmonic, helix, fourier, etc.')

    improper_dihedrals: Optional[Tuple[int, int, int, int, float]] = Field(..., description = 'Dihedral/torsion parameters.')
    improper_dihedrals_pot: Optional[str] = Field(None, description = 'Improper dihedral potential form e.g. harmonic, fourier, etc.')

    charges: Array[float] = Field(..., description = 'Atomic charges.')
    charge_groups: Optional[Array[float]] = Field(None, description = 'Charge groups per atom.')

    nonbonded_pot: Optional[str] = Field(None, description = 'Non-bonded short potential form e.g. Lennard-Jones, Buckingham, etc.')
    nonbonded_params: Array[float] = Field(..., description = 'Non-bonded short potential parameters.')

    exclusions: Optional[str] = Field(None, description = 'Which pairs of bonded atoms to exclude from non-bonded calculations. \
    	The rules to apply in choosing bonded exclusions are specifed in the configuration file using the exclude parameter. The \
    	choices for exclusions are None, 1-2, 1-3, 1-4, etc. With None, no atom pairs are excluded. With 1-2, only atoms that are connected \
    	via a linear bond are excluded. With 1-3, any pair of atoms connected via a bond or bonded to a common third atom are excluded. \
    	With 1-4, any atoms that are connected by a linear bond, or a sequence of two bonds, or a sequence of three bonds, are excluded. \
    	With scaled1-4, exclusions are applied in the same manner as the 1-3 setting, but those pairs that are connected by a sequence of \
    	3 bonds are calculated using the modified 1-4 methods described rather than the standard force calculations.')
    inclusions: Optional[str] = Field(None, description = 'Which pairs of 1-4 bonded atoms to include in non-bonded calculations.')

    name: Optional[str] = Field(None, description = 'Forcefield name e.g. charmm27, amber99, etc.')
    types: Dict[str, str] = Field(None, description = "Atom types specified by the dictionary {'atom': 'type'}.")