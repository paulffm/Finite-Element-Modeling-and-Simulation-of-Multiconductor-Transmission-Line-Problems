# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 11:06:24 2021

.. sectionauthor:: bundschuh
"""

import sys
import gmsh


class GMSHTest:
    """A context manager for gmsh. It initializes and finalizes gmsh.

    Use as in the following example for adding a point:
        >>> with GMSHTest():
        ...    gmsh.model.occ.addPoint(0,0,0)
    """

    def __init__(self, name: str = "test_model"):
        """
        Constructor of the context manager.

        Use as follows:
            >>> with GMSHTest():
                # Your code

        Parameters
        ----------
        name : str, optional
            Name of the model. The default is "test_model".

        Returns
        -------
        None.

        """
        self.name = name

    def __enter__(self):
        gmsh.initialize(sys.argv)

    def __exit__(self, exc_type, exc_value, traceback):
        gmsh.finalize()
