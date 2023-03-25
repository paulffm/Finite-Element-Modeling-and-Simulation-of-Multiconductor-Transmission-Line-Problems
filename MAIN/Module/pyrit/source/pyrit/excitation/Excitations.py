# coding=utf-8
"""
Implementation of Excitations

.. sectionauthor:: Bundschuh
"""
from typing import List, Any
from warnings import warn

from .Exci import Exci
from .FieldCircuitCoupling import FieldCircuitCoupling


class Excitations:
    """Container class that collects all excitations of a problem."""

    def __init__(self, *excis: Exci):
        """
        Constructor of Exciations

        Parameters
        ----------
        excis : Exci
            A number of instances of Exci
        """
        self.__excis = {}
        if len(excis) > 0:
            self.add_exci(*excis)

    def __iter__(self):
        return self.__excis.values().__iter__()

    def add_exci(self, *excis: Exci) -> None:
        """
        Add a number of Excis

        Parameters
        ----------
        excis : Exci
            A number of instances of Exci

        Returns
        -------
        None
        """
        for exci in excis:
            self.__excis[exci.ID] = exci

    def delete_exci(self, *IDs: int) -> None:
        """
        Delete a number of Excis specified by their IDs

        Parameters
        ----------
        IDs : int
            The IDs

        Returns
        -------
        None
        """
        for exci_id in IDs:
            del self.__excis[exci_id]

    def get_exci(self, ID: int) -> Exci:
        """
        Get a specific Exci

        Parameters
        ----------
        ID : int
            The ID of the requested Exci

        Returns
        -------
        exci: Exci
            An Exci
        """
        return self.__excis[ID]

    def get_exci_id_from_name(self, name: str) -> int:
        """
        Get the ID of a specific Exci from its name.

        The ID of the first Exci with the given name is returned. If multiple Exci's with the same name are present,
        a warning is given.

        Parameters
        ----------
        name : str
                The name of the requested Exci

        Returns
        -------
        id: int
            The ID of the Exci with the given name
        """
        ids = self.get_ids()
        exci_names = [self.get_exci(id).name for id in ids]
        if exci_names.count(name) > 1:
            msg = (f"The requested excitation named {name} is listed more than once. Returning the first Exci. Be "
                   f"aware of the global counter ID. Have you added the same excitation more than once or reused its "
                   f"name?")
            warn(msg)
        return ids[exci_names.index(name)]

    def get_exci_from_name(self, name: str) -> Exci:
        """
        Get a specific Exci from its name. The first Exci with the given name is returned.

        Parameters
        ----------
        name : str
                The name of the requested Exci

        Returns
        -------
        exci: Exci
            The Exci with the given name
        """
        return self.get_exci(self.get_exci_id_from_name(name))

    def get_ids(self) -> List[int]:
        """
        Get the IDs of all stored Excis

        Returns
        -------
        ids: List[int]
            A list of all IDs
        """
        return list(self.__excis.keys())

    @property
    def is_linear(self) -> bool:
        """Returns True if all excitations are linear.

        Returns
        -------
        linear : bool
        """
        return all((exci.is_linear for exci in self.__excis.values()))

    @property
    def is_time_dependent(self) -> bool:
        """Returns True if at least one excitation is time dependent.

        Returns
        -------
        time_dependent : bool
        """
        return any((exci.is_time_dependent for exci in self.__excis.values()))

    @property
    def is_homogeneous(self) -> bool:
        """Returns True if all excitations are homogeneous.

        Returns
        -------
        homogeneous : bool
        """
        return all((exci.is_homogeneous for exci in self.__excis.values()))

    @property
    def is_constant(self) -> bool:
        """Returns True if all excitations are constant.

        Returns
        -------
        constant : bool
        """
        return all((exci.is_constant for exci in self.__excis.values()))

    @property
    def is_fcc(self) -> bool:
        """Returns True if it contains any exci of type FieldCircuitCoupling"""
        return any((isinstance(exci, FieldCircuitCoupling) for exci in self.__excis.values()))

    def update(self, name: str, val: Any):
        """Update a value to a name.

        Calls the function :py:meth:`pyrit.excitation.Exci.update` on every excitation.

        Parameters
        ----------
        name : str
            The name of the field to update.
        val : Any
            The value to update.
        """
        for exci in self.__excis.values():
            exci.update(name, val)
