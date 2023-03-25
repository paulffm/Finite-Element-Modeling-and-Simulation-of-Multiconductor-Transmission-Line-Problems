# coding=utf-8
"""Charge density.

.. sectionauthor:: Bundschuh
"""

from . import Exci


class ChargeDensity(Exci):
    """Class representing a charge density."""

    def __init__(self, value, name: str = ''):
        super().__init__(value, name)
