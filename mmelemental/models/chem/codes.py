from mmelemental.models.base import Base
from pydantic import validator, Field, ValidationError
from mmelemental.util.decorators import req_rdkit
import os

try:
    from rdkit import rdBase, Chem
    from rdkit.Chem import AllChem
    if not os.environ.get('debugMMC'):
        rdBase.DisableLog('rdApp.error')
except ImportError:
    Chem = AllChem = None
        
class ChemCode(Base):
    code: str = Field(
        ...,
        description = 'A chemical code that describes a molecule or molecular pattern e.g. smiles, smarts, etc. '
        'See :class: ``Identifiers`` for supported codes.'
    )

    class _CodesSupported:
        codes = ('Smiles', 'Smarts', 'Inchi', 'FASTA', 'HELM', 'Sequence')

    @req_rdkit
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @req_rdkit
    @validator('code')
    def valid_code(cls, code):
        for ccode in ChemCode._CodesSupported.codes:
            function = getattr(Chem, f"MolFrom{ccode}")
            if function(code):
                return code
        raise ValidationError

    @property
    def code_type(self):
        for ccode in ChemCode._CodesSupported.codes:
            function = getattr(Chem, f"MolFrom{ccode}")
            if function(self.code):
                return ccode
        raise ValidationError