from pydantic import Field, constr, validator
from mmelemental.models.base import ProtoModel
import qcelemental
from .di_params import DihedralsParams
from typing import Optional, Dict, Any

__all__ = ["Dihedrals"]


class Dihedrals(ProtoModel):
    params: DihedralsParams = Field(..., description="Dihedral parameters model.")
    angles: qcelemental.models.types.Array[float] = Field(
        None, description="Equilibrium dihedral angles. Default unit is degrees."
    )
    angles_units: Optional[str] = Field(
        "degrees", description="Equilibrium dihedral angle units."
    )

    # Constructors
    @classmethod
    def from_file(
        cls,
        filename: str,
        dtype: Optional[str] = None,
        translator: Optional[str] = None,
        **kwargs,
    ) -> "Dihedrals":
        """
        Constructs a Dihedrals object from a file.
        Parameters
        ----------
        filename: str
            The topology or FF filename to build from.
        dtype: Optional[str], optional
            The type of file to interpret e.g. psf. If unset, mmelemental attempts to discover the file type.
        translator: Optional[str], optional
            Translator name e.g. mmic_parmed. Takes precedence over dtype. If unset, MMElemental attempts 
            to find an appropriate translator if it is registered in the :class:``TransComponent`` class. 
        **kwargs: Optional[Dict[str, Any]], optional
            Any additional keywords to pass to the constructor.
        Returns
        -------
        Dihedrals
            A constructed Dihedrals object.
        """

        fileobj = FileOutput(path=filename)
        dtype = dtype or fileobj.ext.strip(".")
        ext = "." + dtype

        if translator:
            raise NotImplementedError

        with open(filename, "r") as fp:
            data = json.load(fp)

        return cls(**data)

    @classmethod
    def from_data(cls, data: Any, **kwargs) -> "Dihedrals":
        """
        Constructs a Dihedrals object from a data object.
        Parameters
        ----------
        data: Any
            Data to construct Dihedrals from.
        **kwargs: Optional[Dict[str, Any]], optional
            Additional kwargs to pass to the constructors.
        Returns
        -------
        Dihedrals
            A constructed Dihedrals object.
        """
        if hasattr(data, "to_schema"):
            return data.to_schema(**kwargs)
        else:
            raise NotImplementedError

    def to_file(
        self,
        filename: str,
        dtype: Optional[str] = None,
        translator: Optional[str] = None,
        **kwargs: Dict[str, Any],
    ) -> None:
        """Writes the Dihedrals to a file.
        Parameters
        ----------
        filename : str
            The filename to write to
        dtype : Optional[str], optional
            The type of file to write (e.g. psf, top, etc.), attempts to infer dtype from
            file extension if not provided.
        translator: Optional[str], optional
            Translator name e.g. mmic_parmed. Takes precedence over dtype. If unset, MMElemental attempts 
            to find an appropriate translator if it is registered in the :class:``TransComponent`` class. 
        **kwargs: Optional[str, Dict], optional
            Additional kwargs to pass to the constructor.
        """
        if not dtype:
            from pathlib import Path

            ext = Path(filename).suffix
        else:
            ext = "." + dtype

        if translator or ext != ".json":
            raise NotImplementedError

        stringified = self.json(**kwargs)
        mode = kwargs.pop("mode", "w")

        with open(filename, mode) as fp:
            fp.write(stringified)

    def to_data(
        self,
        dtype: Optional[str] = None,
        translator: Optional[str] = None,
        **kwargs: Dict[str, Any],
    ) -> "ToolkitModel":
        """
        Constructs a toolkit-specific data object from MMSchema Dihedrals.
        Which toolkit-specific component is called depends on which package is installed on the system.
        Parameters
        ----------
        translator: Optional[str], optional
            Translator name e.g. mmic_parmed. Takes precedence over dtype. If unset, MMElemental attempts 
            to find an appropriate translator if it is registered in the :class:``TransComponent`` class. 
        dtype: Optional[str], optional
            Data type e.g. MDAnalysis, parmed, etc.
        **kwargs: Optional[Dict[str, Any]]
            Additional kwargs to pass to the constructors.
        Results
        -------
        ToolkitModel
            Toolkit-specific Dihedrals object
        """
        raise NotImplementedError

    def __eq__(self, other):
        """
        Checks if two molecules are identical. This is a molecular identity defined
        by scientific terms, and not programing terms, so it's less rigorous than
        a programmatic equality or a memory equivalent `is`.
        """

        if isinstance(other, dict):
            other = Dihedrals(**other)
        elif isinstance(other, Dihedrals):
            pass
        else:
            raise TypeError(f"Comparison not understood of type '{type(other)}'.")

        return self.get_hash() == other.get_hash()

    # Validators
    @validator("angles")
    def _angles_length(cls, v, values):
        assert len(v.shape) == 1, "Dihedral angles must be a 1D array!"
        return v

    # Propreties
    @property
    def hash_fields(self):
        return ["charges", "params"]

    @property
    def form(self):
        return self.params.__class__.__name__

    def get_hash(self):
        """
        Returns the hash of the force field object.
        """

        m = hashlib.sha1()
        concat = ""

        # np.set_printoptions(precision=16)
        for field in self.hash_fields:
            data = getattr(self, field)
            if data is not None and field != "params":
                # if field == "nonbonded":
                #    data = qcelemental.models.molecule.float_prep(data, GEOMETRY_NOISE)

                concat += json.dumps(data, default=lambda x: x.ravel().tolist())

        m.update(concat.encode("utf-8"))
        return m.hexdigest()

    def dict(self, *args, **kwargs):
        kwargs["exclude"] = {"provenance"}
        return super().dict(*args, **kwargs)
