""" A list of helpful decorators for MMElemental """

__all__ = ["deprecated", "require"]

import warnings
import functools
import importlib


def require(tk_name):
    def inner_require(func):
        @functools.wraps(func)
        def inner_func(*args, **kwargs):
            if not importlib.import_module(tk_name):
                raise ModuleNotFoundError(f"Could not find or import {tk_name}.")
            return func(*args, **kwargs)

        return inner_func

    return inner_require


def deprecated(func):
    @functools.wraps(func)
    def func_inner(*args, **kwargs):
        warnings.warn(
            f"{func.__name__} is deprecated and will be removed in the future",
            DeprecationWarning,
        )
        return func(*args, **kwargs)

    return func_inner
