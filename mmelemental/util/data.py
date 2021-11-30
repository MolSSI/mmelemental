import numpy
from typing import Any, Union

# Numpy dtypes for heterogenous arrays
NUMPY_INT = "i8"
NUMPY_FLOAT = "f8"
NUMPY_UNI = "U4"

# Rounding quantities for hashing
GEOMETRY_NOISE = 8
VELOCITY_NOISE = 8
FORCE_NOISE = 8
MASS_NOISE = 6
CHARGE_NOISE = 4


def float_prep(array: numpy.ndarray, around: int) -> numpy.ndarray:
    """
    Rounds floats to a common value and build positive zeros to prevent hash conflicts.
    """
    if isinstance(array, (list, numpy.ndarray)):
        # Round array
        array = numpy.around(array, around)
        # Flip zeros
        array[numpy.abs(array) < 5 ** (-(around + 1))] = 0

    elif isinstance(array, (float, int)):
        array = round(array, around)
        if array == -0.0:
            array = 0.0
    else:
        raise TypeError("Type '{}' not recognized".format(type(array).__name__))

    return array


def get_dtype(data: Any) -> Union[str, None]:
    """Returns a numpy data type for an object. Integers and floats are represented in 8 bytes (64 bits)
    precision while strings in 4 characters (128 bits) UTF-8 unicode.

    Returns
    -------
    str or None
        Returns data type if found otherwise None.

    Raises
    ------
    ValueError
        When `data` type is not a str, float, int, list, tuple, or numpy.ndarray.
    """
    if isinstance(data, str):
        return NUMPY_UNI
    elif isinstance(data, float):
        return NUMPY_FLOAT
    elif isinstance(data, int):  # check for more types?
        return NUMPY_INT
    elif isinstance(data, (list, tuple, numpy.ndarray)):
        if isinstance(data[0], (list, tuple, numpy.ndarray)):
            if isinstance(data[0], numpy.ndarray):
                if data.dtype != object:
                    return None  # homogenous array
            dtype = ",".join([get_dtype(item) for item in data[0]])
        elif all([isinstance(item, str) for item in data]):
            dtype = str
        else:
            return None  # We assume homogenous array in MMSchema esle dtype=object
        return dtype
    else:
        raise ValueError(f"Data type {type(data)} not supported.")
