# coding=utf-8
"""Current density

.. sectionauthor: Bundschuh
"""

from . import Exci


class CurrentDensity(Exci):
    """Class representing a current density."""

    def __init__(self, value, name: str = ''):
        super().__init__(value, name)
