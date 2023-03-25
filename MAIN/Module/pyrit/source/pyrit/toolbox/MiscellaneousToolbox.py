# coding=utf-8
"""
Toolbox for miscellaneous functions

.. sectionauthor:: bundschuh
"""

from inspect import signature
from typing import Callable, Any, List
from time import time
from contextlib import contextmanager
import cProfile
import pstats


def evaluate_function(fun: Callable[..., Any], dependencies: List[str] = None, **kwargs) -> Any:
    """
    Evaluates the function fun. The needed parameters are taken from kwargs.

    Parameters
    ----------
    fun : Callable[..., Any]
        The function to evaluate
    dependencies : List[str], optional
        A list of dependencies. If not given it is calculated internally.
    kwargs :
        All Attributes

    Returns
    -------
    The return value of the function fun
    """
    if dependencies is None:
        dependencies = list(signature(fun).parameters.keys())
    kw = {}
    for key in dependencies:
        if key in kwargs:
            kw[key] = kwargs[key]
        else:
            raise Exception("Key not in dependencies")

    return fun(**kw)


def timer_func(func):
    """Decorator that measures the time for a function."""

    def wrap_func(*args, **kwargs):
        t1 = time()
        result = func(*args, **kwargs)
        t2 = time()
        print(f'Function {func.__name__!r} executed in {((t2 - t1) * 1000):.9f} ms')
        return result

    return wrap_func


@contextmanager
def profile(sort_key=None):
    """A profiler context manager.

    Collects the execution times of functions.

    Parameters
    ----------
    sort_key: str
        The key to sort the results. Eiter 'time', 'cumulative', 'calls', 'module' or any other `SortKey` from `pstats`
        module.
    """
    sort_keys = {'time': pstats.SortKey.TIME, 'cumulative': pstats.SortKey.CUMULATIVE, 'calls': pstats.SortKey.CALLS,
                 'module': pstats.SortKey.FILENAME}
    if sort_key is None:
        sort_key = 'time'
    sort_key = sort_keys.setdefault(sort_key, sort_key)

    with cProfile.Profile() as pr:
        yield pr
    stats = pstats.Stats(pr)
    stats.sort_stats(sort_key)
    stats.print_stats()
