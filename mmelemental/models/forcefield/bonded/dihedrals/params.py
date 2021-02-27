from pydantic import root_validator
from mmelemental.models.base import ProtoModel

__all__ = ["BondsParams"]

supported_potentials = ("Harmonic",)


class BondsParams(ProtoModel):
    @root_validator
    def _is_registered(cls, values):
        if cls.__name__ not in supported_potentials:
            raise NotImplementedError(
                f"{cls.__name__} is not supported in MMElemental."
            )
        return values

    def dict(self, *args, **kwargs):
        kwargs["exclude"] = {"provenance"}
        return super().dict(*args, **kwargs)
