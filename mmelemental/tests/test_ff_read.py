from mmelemental.models import forcefield as ff
import glob
import mm_data

json_files = [mm_data.ffs["simple.json"]]


def pytest_generate_tests(metafunc):
    if "json_file" in metafunc.fixturenames:
        metafunc.parametrize("json_file", json_files)


def test_read_json(json_file):
    ff.ForceField.from_file(json_file)
