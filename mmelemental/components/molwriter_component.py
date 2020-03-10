try:
    from rdkit import rdBase, Chem
    from rdkit.Chem import AllChem
except:
    raise ModuleNotFoundError('Make sure rdkit is installed for code validation.')

from mmcomponents.components.blueprints.generic_component import GenericComponent
from mmelemental.models.utils.output import FileOutput

class MMoleculeWriter(GenericComponent):

    @classmethod
    def input(cls):
        return MMolecule

    @classmethod
    def output(cls):
        return FileOutput
