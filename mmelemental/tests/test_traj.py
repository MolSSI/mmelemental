import pytest
from mmelemental.models import Trajectory, Molecule
from cmselemental.util import yaml_import, which_import
import numpy
from pathlib import Path
import mm_data
import os

using_yaml = pytest.mark.skipif(
    yaml_import() is False,
    reason="Not detecting module pyyaml or ruamel.yaml. Install package if necessary and add to envvar PYTHONPATH",
)

using_nglview = pytest.mark.skipif(
    which_import("nglview", return_bool=True) is False,
    reason="Not detecting module nglview. Install package if necessary and add to envvar PYTHONPATH",
)

file_extensions = [
    ".json",
    pytest.param(".yaml", marks=using_yaml),
]


@pytest.mark.parametrize(
    "data",
    [
        {
            "geometry": numpy.random.rand(10 * 3 * 2),
            "natoms": 10,
            "ndim": 3,
            "nframes": 2,
            "timestep": 1,
        },
        {
            "natoms": 10,
            "ndim": 3,
            "nframes": 4,
            "timestep": 1,
        },
        {
            "geometry": numpy.random.rand(20 * 2 * 3),
            "velocities": numpy.random.rand(20 * 2 * 3),
            "forces": numpy.random.rand(20 * 2 * 3),
            "natoms": 20,
            "ndim": 2,
            "nframes": 3,
            "timestep": 0.1,
        },
    ],
)
def test_mmelemental_traj(data):
    traj = Trajectory(**data)
    fpath = Path("tmp.json")
    traj.to_file(fpath.name, indent=4)
    assert fpath.is_file()
    fpath.unlink()


@pytest.mark.parametrize(
    "mfile",
    [mol for mol in mm_data.mols.values() if mol.endswith(".json")],
)
def test_mmelemental_toptraj(mfile):
    mol = Molecule.from_file(mfile)
    top = mol.get_topology()
    natoms = len(top.symbols)
    nframes = 20
    ndim = 3

    data = {
        "top": top,
        "geometry": numpy.random.rand(natoms * ndim * nframes),
        "nframes": nframes,
        "ndim": ndim,
        "natoms": natoms,
        "timestep": 0.2,
    }
    traj = Trajectory(**data)
    fpath = Path("tmp.json")
    traj.to_file(fpath.name, indent=4)
    assert fpath.is_file()
    fpath.unlink()


@pytest.mark.parametrize("ext", file_extensions)
def test_mmelemental_traj_files(ext):
    for filen, filep in mm_data.trajs.items():
        if ext in filen:
            mm_traj = Trajectory.from_file(filep, all_frames=True)
            assert isinstance(mm_traj, Trajectory)

            fpath = Path(f"tmp{ext}")
            mm_traj.to_file(fpath.name)
            assert fpath.is_file()
            fpath.unlink()
