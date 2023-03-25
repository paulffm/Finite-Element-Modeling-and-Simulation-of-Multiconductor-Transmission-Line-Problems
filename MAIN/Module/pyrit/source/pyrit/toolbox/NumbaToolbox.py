# coding=utf-8
"""Toolbox for numba"""
import os
from numba import njit


def njit_cache(*args, **kwargs):
    """Sets the cach keyword of numba njit decorator depending on the environment variable "NUMBA_DISABLE_JIT"."""
    # pylint: disable=invalid-name
    NUMBA_DISABLE_JIT = os.environ.get("NUMBA_DISABLE_JIT")
    if NUMBA_DISABLE_JIT is not None:
        if NUMBA_DISABLE_JIT == '1':
            return njit(*args, cache=False, **kwargs)
    return njit(*args, cache=True, **kwargs)


def diagnose(func):
    """Decorator that calls diagnostics.

    This decorator runs diagnostics for functions with the jit decorator from numba with the keyword `parallel=True`.
    """
    def wrapper(*args, **kwargs):
        val = func(*args, **kwargs)
        func.parallel_diagnostics()
        func.inspect_types()
        return val

    return wrapper
