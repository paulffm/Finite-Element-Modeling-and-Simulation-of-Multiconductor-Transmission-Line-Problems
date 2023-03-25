# coding=utf-8
"""A specialised logger for `Pyrit`.

The :py:mod:`logging` module from the python standard library is used for logging. For general information, see that
documentation. The colored log-messages are realized with the module :py:mod:`colorama`. For each type of log-message
there is an attribute that determines the used color.
This module provides functions to change single colors, the level for logging and the format string.

.. note::

    If you want to change a default value, this has to be done **bevore** other modules from `Pyrit` are loaded.

Example
-------
Example of changing the default color for log-messages of type INFO

>>> import pyrit
>>> pyrit.set_color('INFO',colorama.Fore.WHITE)
>>> from pyrit.geometry import Geometry

How to use this module
----------------------
It is intended to use this module as follows. At the top of a file, a new logger is created. This logger is then used
for this whole file with the standard methods for logging:

>>> from pyrit import get_logger
>>> logger = get_logger(__name__)
>>> ...
>>> logger.info("This is an info log-message")
>>> logger.warning("This is a warning")

.. sectionauthor:: Bundschuh
"""
# pylint: disable=invalid-name, global-statement

import logging
from typing import Any
import colorama

colorama.init(autoreset=True)

LEVEL = logging.WARNING  #: The level for the logging module

COLOR_DEBUG = colorama.Fore.WHITE  #: Color for debug log-messages
COLOR_INFO = colorama.Fore.GREEN  #: Color for info log-messages
COLOR_WARNING = colorama.Fore.YELLOW  #: Color for warning log-messages
COLOR_ERROR = colorama.Fore.RED  #: Color for error log-messages
COLOR_CRITICAL = colorama.Fore.RED + colorama.Style.DIM  #: Color for critical log-messages

#: The format string for log-messages
LOG_FORMAT = '%(levelname)s: Module %(module)s (Line %(lineno)s) :: %(message)s'

FORMATS = {}


def reset_colors():
    """Reset all colors to the default values"""
    default_formats = {
        logging.DEBUG: COLOR_DEBUG + LOG_FORMAT,
        logging.INFO: COLOR_INFO + LOG_FORMAT,
        logging.WARNING: COLOR_WARNING + LOG_FORMAT,
        logging.ERROR: COLOR_ERROR + LOG_FORMAT,
        logging.CRITICAL: COLOR_CRITICAL + LOG_FORMAT
    }
    global FORMATS
    FORMATS = default_formats.copy()


reset_colors()


def set_level(level: Any):
    """Set the level for the logging output.

    Parameters
    ----------
    level : Any
        The level. See :py:mod:`logging`.
    """
    global LEVEL
    LEVEL = level


def set_color(level: Any, color: str):
    """Set a color for the given level.

    If the level is not known, it is ignored.

    Parameters
    ----------
    level : Any
        The level. See :py:mod:`logging`.
    color : str
        The color for the level. See the module :py:mod:`colorama`.
    """
    # pylint: disable=global-variable-not-assigned
    global FORMATS
    try:
        FORMATS[level] = color + LOG_FORMAT
    except KeyError:
        try:
            FORMATS[logging.getLevelName(level)] = color + LOG_FORMAT
        except KeyError:
            print(f"The level '{level}' is not known. It is ignored.")


def disable_colors():
    """Disable colored output"""
    keys = [logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL]
    for key in keys:
        set_color(key, '')


def set_format(new_format: str):
    """Set the format string for log-messages.

    Parameters
    ----------
    new_format : str
        The new format string. See :py:mod:`logging` for for information.
    """
    global LOG_FORMAT
    LOG_FORMAT = new_format


class CustomFormatter(logging.Formatter):
    """A custom formatter that uses the colored output."""

    def format(self, record):
        log_fmt = FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)


def get_logger(name: str) -> logging.Logger:
    """Get a logger object

    Parameters
    ----------
    name : str
        The name of the logger. It is recommended to use __name__.

    Returns
    -------
    logger : logging.Logger
        The logger object.
    """
    logger = logging.getLogger(name)
    ch = logging.StreamHandler()
    ch.setFormatter(CustomFormatter())
    # ch.setLevel('INFO')
    logger.setLevel('INFO')
    logger.addHandler(ch)

    return logger
