from mmelemental.models.base import ProtoModel

__all__ = ["ToolkitMol"]


class ToolkitMol(ProtoModel):
    """A generic base class that acts as a wrapper for toolkit molecules
    TODO: Delete this class and move check_name() to wherever it's needed (rdkit?!)"""

    class Config:
        arbitrary_types_allowed = True

    @staticmethod
    def check_name(name) -> str:
        """Returns atom name of langth 4 characters."""
        assert len(name) <= 4

        if len(name) != 4:
            if len(name) == 1:
                name = " " + name + "  "
            elif len(name) == 2:
                name = " " + name + " "
            elif len(name) == 3:
                name = " " + name
        return name
