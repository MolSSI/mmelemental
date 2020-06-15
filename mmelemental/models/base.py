from qcelemental import models
from pydantic import Field

from typing import Dict
from qcelemental.extras import get_information
from qcelemental.models.common_models import Provenance

def provenance_stamp(routine: str) -> Dict[str, str]:
    """Return dictionary satisfying QCSchema,
    https://github.com/MolSSI/QCSchema/blob/master/qcschema/dev/definitions.py#L23-L41
    with QCElemental's credentials for creator and version. The
    generating routine's name is passed in through `routine`.
    """
    return {"creator": "MMElemental", "version": get_information("version"), "routine": routine}

class Base(models.ProtoModel):
    provenance: Provenance = Field(
        provenance_stamp(__name__),
        description="The provenance information about how this object (and its attributes) were generated, "
        "provided, and manipulated.",
    )