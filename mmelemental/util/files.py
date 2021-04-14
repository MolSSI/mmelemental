import uuid


def random_file(suffix=""):
    return str(uuid.uuid4()) + suffix
