""" A list of helpful decorators for MMElemental """

__all__ = ["req_openmm", "require", "req_mda", "req_parmed"]

import warnings
import functools
import importlib

try:
    import simtk.openmm as mm
    import simtk.openmm.app as app

    HAS_OPENMM = True
    try:
        from simtk.openmm.app.internal import unitcell
    except ImportError:
        unitcell = None
        SUPPORTED_VERSION = False
    else:
        SUPPORTED_VERSION = True
except ImportError:
    HAS_OPENMM = False
else:
    del mm, app, unitcell

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

try:
    import MDAnalysis

    HAS_MDA = True
except ImportError:
    HAS_MDA = False

try:
    import parmed

    HAS_PARMED = True
except ImportError:
    HAS_PARMED = False


def req_parmed(func):
    global HAS_PARMED

    @functools.wraps(func)
    def func_inner(*args, **kwargs):
        if not HAS_PARMED:
            raise ModuleNotFoundError("Could not find or import parmed.")
        return func(*args, **kwargs)

    return func_inner


def req_mda(func):
    global HAS_MDA

    @functools.wraps(func)
    def func_inner(*args, **kwargs):
        if not HAS_MDA:
            raise ModuleNotFoundError("Could not find or import MDAnalysis.")
        return func(*args, **kwargs)

    return func_inner


def require(tk_name):
    def inner_require(func):
        @functools.wraps(func)
        def inner_func(*args, **kwargs):
            if not importlib.import_module(tk_name):
                raise ModuleNotFoundError(f"Could not find or import {tk_name}.")
            return func(*args, **kwargs)

        return inner_func

    return inner_require


def req_openmm(func):
    global HAS_OPENMM

    @functools.wraps(func)
    def func_inner(*args, **kwargs):
        if not HAS_OPENMM:
            raise ModuleNotFoundError("Could not find or import OpenMM.")
        if not SUPPORTED_VERSION:
            raise ImportError("You must have at least OpenMM 6.3 installed.")
        return func(*args, **kwargs)

    return func_inner


def deprecated(func):
    @wraps(func)
    def func_inner(*args, **kwargs):
        warnings.warn(
            f"{func.__name__} is deprecated and will be removed in the future",
            DeprecationWarning,
        )
        return func(*args, **kwargs)

    return func_inner
