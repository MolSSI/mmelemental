from mmelemental.models import forcefield as ff
import glob
import os
from mm_data import data_ff_dir

json_files = glob.glob(os.path.join(data_ff_dir, "*.json"))


def pytest_generate_tests(metafunc):
    if "json_file" in metafunc.fixturenames:
        metafunc.parametrize("json_file", json_files)


def test_read_json(json_file):
    ff.ForceField.from_file(json_file)
