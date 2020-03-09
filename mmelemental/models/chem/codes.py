from qcelemental import models
from pydantic import validator, Field, ValidationError
import os

try:
    from rdkit import rdBase, Chem
    from rdkit.Chem import AllChem
    if not os.environ.get('debugMMC'):
        rdBase.DisableLog('rdApp.error')
    rdkAvail = True
except Exception:
    rdkAvail = False
        
class ChemCode(models.ProtoModel):
    code: str

    class _CodesSupported:
        codes = ('Smiles', 'Smarts', 'Inchi', 'FASTA', 'HELM', 'Sequence')

    @validator('code')
    def validCode(cls, code):
        if rdkAvail:
            for ccode in ChemCode._CodesSupported.codes:
                function = getattr(Chem, f"MolFrom{ccode}")
                if function(code):
                    return code
            raise ValidationError
        else:
            raise ModuleNotFoundError('Make sure rdkit is installed for code validation.')

    @property
    def codeType(self):
        if rdkAvail:
            for ccode in ChemCode._CodesSupported.codes:
                function = getattr(Chem, f"MolFrom{ccode}")
                if function(self.code):
                    return ccode
            raise ValidationError
        else:
            raise ModuleNotFoundError('Make sure rdkit is installed for code validation.')