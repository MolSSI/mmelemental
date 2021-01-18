from mmic.components.blueprints.generic_component import GenericComponent
from mmelemental.models.molecule.parmed_molecule import ParmedMolecule
from mmelemental.models.molecule.io_molecule import MolInput
from mmelemental.models.molecule.mm_molecule import Molecule
from mmelemental.models.molecule.gen_molecule import ToolkitMolecule
from mmelemental.models.forcefield.mm_ff import ForceField
from mmelemental.models.forcefield.io_ff import FFInput
from typing import Dict, Any, List, Tuple, Optional

try:
    from parmed import structure, topologyobjects
except:
    raise ModuleNotFoundError('Make sure parmed is installed.')

class MolToParmedComponent(GenericComponent):
    """ A component for converting Molecule to ParmEd molecule object. """
    @classmethod
    def input(cls):
        return Molecule

    @classmethod
    def output(cls):
        return ParmedMolecule

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, ParmedMolecule]:

        mmol = inputs
        pmol = structure.Structure()

        for index, symb in enumerate(mmol.symbols):

            name = mmol.names[index]
            #name = ToolkitMolecule.check_name(name)

            atomic_number = mmol.atomic_numbers[index]
            mass = mmol.masses[index]

            # Will likely lose FF-related info ... but then Molecule is not supposed to store any params specific to FFs
            atom = topologyobjects.Atom(list=None, atomic_number=atomic_number, name=name, type=symb,
                 mass=mass, nb_idx=0, solvent_radius=0.0,
                 screen=0.0, tree='BLA', join=0.0, irotat=0.0, occupancy=1.0,
                 bfactor=0.0, altloc='', number=-1, rmin=None, epsilon=None,
                 rmin14=None, epsilon14=None)

            resname, resnum = mmol.residues[index]
            #classparmed.Residue(name, number=- 1, chain='', insertion_code='', segid='', list=None)[source]

            pmol.add_atom(atom, resname, resnum + 1, chain='', inscode='', segid='')

        for i, j , btype in mmol.connectivity:
            pmol.atoms[i].bond_to(pmol.atoms[j])

        pmol.coordinates = mmol.geometry

        return True, ParmedMolecule(mol=pmol)

class ParmedToMolComponent(GenericComponent):
    """ A component for converting ParmEd molecule to Molecule object. """

    @classmethod
    def input(cls):
        return MolInput

    @classmethod
    def output(cls):
        return Molecule

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None) -> Tuple[bool, Dict[str, Any]]:
        
        if inputs.data:
            dtype = inputs.data.dtype
            assert dtype == 'parmed'
            pmol = inputs.data
        elif inputs.code:
            raise NotImplementedError('ParmEd does not support instantiating molecule objects from chemical codes.')          
        elif inputs.file:
            dtype = '.' + inputs.file.ext
            from mmelemental.models.molecule.parmed_molecule import ParmedMolecule
            pmol = ParmedMolecule.build(inputs, dtype)
        else:
            raise NotImplementedError(f'Data type {dtype} not yet supported.')
        
        if inputs.args:
            orient = inputs.args.get('orient')
            validate = inputs.args.get('validate')
            kwargs = inputs.args.get('kwargs')
        else:
            orient, validate, kwargs = False, None, None

        symbs = [atom.element_name for atom in pmol.mol.atoms]
        names = [atom.name for atom in pmol.mol.atoms]

        if hasattr(pmol.mol.atoms[0], 'mass'):
            masses = [atom.mass for atom in pmol.mol.atoms]

        if hasattr(pmol.mol.atoms[0], 'residue'):
            residues = [(atom.residue.name, atom.residue.idx) for atom in pmol.mol.atoms]

        connectivity = []

        for bond in pmol.mol.bonds:
            connectivity.append((bond.atom1.idx, bond.atom2.idx, bond.order))

        geo = pmol.mol.coordinates

        input_dict = {'symbols': symbs, 
                      'geometry': geo, 
                      'residues': residues, 
                      'connectivity': connectivity,
                      'masses': masses,
                      'names': names}        

        if kwargs:
            input_dict.update(kwargs)

        if inputs.code:
            return True, Molecule(orient=orient, validate=validate, identifiers={dtype: inputs.code}, **input_dict)
        else:
            return True, Molecule(orient=orient, validate=validate, **input_dict)

class ParmedToFFComponent(GenericComponent):
    """ A component for converting ParmEd top object to MMSchema ForceField (FF) object. """

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
        timeout: Optional[int] = None) -> Tuple[bool, Dict[str, Any]]:
        
        if inputs.data:
            dtype = inputs.data.dtype
            assert dtype == 'parmed'
            ff = inputs.data
        elif inputs.file:
            dtype = '.' + inputs.file.ext
            from mmelemental.models.forcefield.parmed_ff import ParmedFF
            ff = ParmedFF.build(inputs, dtype)
        else:
            raise NotImplementedError(f'Data type {dtype} not yet supported.')

        if inputs.args:
            if 'software' in inputs.args:
                if inputs.args['software'] == 'gromacs':
                    if len(ff.ff.itps) > 1:
                        raise ValueError(f'More than 1 forcefield name detected: {ff.ff.itp}.')
                    else:
                        print(f"Forcefield {ff.ff.itps[0].split('/')[0].split('.ff')[0]} detected.")

        data = [(atom.name, atom.charge) for atom in ff.ff.atoms]
        names, charges = zip(*data)
        params = [(atom.rmin, atom.rmin_14, atom.epsilon, atom.epsilon_14) for atom in ff.ff.atoms]

        bonds = [(bond.type.req, bond.type.k) for bond in ff.ff.bonds]
        bonds_type = [bond.funct for bond in ff.ff.bonds]

        angles = [(angle.type.theteq, angle.type.k) for angle in ff.ff.angles]
        angles_type = [angle.funct for angle in ff.ff.angles]

        # Why the hell does parmed interpret dihedral.type as a list sometimes? See dialanine.top, last 8 dihedral entries. WEIRD!
        dihedrals = [(dihedral.type[0].phase, dihedral.type[0].phi_k, dihedral.type[0].per, \
                    dihedral.type[0].scee, dihedral.type[0].scnb) if isinstance(dihedral.type, list) \
                    else (dihedral.type.phase, dihedral.type.phi_k, dihedral.type.per, \
                    dihedral.type.scee, dihedral.type.scnb) for dihedral in ff.ff.dihedrals]
        dihedrals_type = [dihedral.funct for dihedral in ff.ff.dihedrals]

        improper_dihedrals = None
        charge_groups = None
        exclusions = ff.ff.nrexcl
        inclusions = None

        input_dict = {'bonds': bonds,
                      'angles': angles,
                      'dihedrals': dihedrals,
                      'improper_dihedrals': improper_dihedrals,
                      'charges': charges,
                      'charge_groups': charge_groups,
                      'nonbonded': params,
                      'exclusions': exclusions,
                      'inclusions': inclusions,
                      'types': names, 
                      'name': 'amber99',
                    }

        return True, ForceField(**input_dict)