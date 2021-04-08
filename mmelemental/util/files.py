import tempfile

def random_file(suffix=""):
    fp = tempfile.NamedTemporaryFile(suffix=suffix)
    fp.close()
    return fp.name