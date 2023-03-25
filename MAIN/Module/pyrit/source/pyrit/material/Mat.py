# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 09:44:38 2021

Contains the class Mat

.. sectionauthor:: bundschuh
"""
# pylint: disable=unused-private-member
from typing import Type, Union, Any

from pyrit import get_logger
from . import MatProperty, Conductivity, Density, DifferentialConductivity, DifferentialPermittivity, Permeability, \
    Permittivity, Reluctivity, Resistivity, ThermalConductivity, VolumetricHeatCapacity

logger = get_logger(__name__)


class PropertyConversionException(Exception):
    """Exception when it comes to an error during the conversion of MatProperty objects."""


class PropertyConverter:
    """Class for converting MatProperties."""

    all_property_classes = [Conductivity, Density, DifferentialConductivity, DifferentialPermittivity, Permeability,
                            Permittivity, Reluctivity, Resistivity, ThermalConductivity, VolumetricHeatCapacity]

    def __init__(self):
        self.convert_functions = {Conductivity: self.from_conductivity,
                                  Density: self.from_density,
                                  DifferentialConductivity: self.from_differential_conductivity,
                                  DifferentialPermittivity: self.from_differential_permittivity,
                                  Permeability: self.from_permeability,
                                  Permittivity: self.from_permittivity,
                                  Reluctivity: self.from_reluctivity,
                                  Resistivity: self.from_resistivity,
                                  ThermalConductivity: self.from_thermal_conductivity,
                                  VolumetricHeatCapacity: self.from_volumetric_heat_capacity}
        self.property_relations = self.compute_property_relations()

    def convertible(self, from_object: MatProperty, to_class) -> bool:
        """True if the `from_object` is convertible to the `to_class` ."""
        return to_class in self.property_relations[type(from_object)]

    def compute_property_relations(self):
        """Computes all possible converting relations

        Returns
        -------
        property_relations : Dict
            Key is a subclass of MatProperty and
            value is a list of subclasses of MatProperty the key can be converted to
        """
        property_relations = {}
        for prop_class_prop in self.all_property_classes:
            prop = prop_class_prop(1, 'dummy')
            classes = []
            for prop_class in self.all_property_classes:
                if prop_class is prop_class_prop:
                    continue
                try:
                    self._convert(prop, prop_class)
                    classes.append(prop_class)
                except PropertyConversionException:
                    pass
            property_relations[prop_class_prop] = classes
        del prop
        return property_relations

    def _convert(self, from_object: MatProperty, to_class) -> MatProperty:
        """Internal convert method without any checks.

        Parameters
        ----------
        from_object : MatProperty
            An instance of a subclass of MatProperty
        to_class :
            A subclass of MatProperty the `from_object` shall be converted to.

        Returns
        -------
        new_prop : MatProperty
            An instance of the `to_class` if the conversion was possible and None if not
        """
        try:
            return self.convert_functions[type(from_object)](from_object, to_class)
        except KeyError as e:
            logger.error("The class ob the MatProperty object %s is not known. Can't convert.", type(from_object))
            raise e

    def convert(self, from_object: MatProperty, to_class) -> MatProperty:
        """Convert a MatProperty to another MatProperty

        Parameters
        ----------
        from_object : MatProperty
            An instance of a subclass of MatProperty
        to_class :
            A subclass of MatProperty the `from_object` shall be converted to.

        Returns
        -------
        new_prop : MatProperty
            An instance of the `to_class` if the conversion was possible

        Raises
        ------
        PropertyConversionException: If conversion is not possible
        """
        return self._convert(from_object, to_class)

    @staticmethod
    def from_conductivity(prop: Conductivity, prop_class):
        """Convert a MatProperty from a conductivity

        Parameters
        ----------
        prop : Conductivity
            Existing MatProperty.
        prop_class :
            Desired output property

        Returns
        -------
        new_prop :
            An object of the desired MatProperty

        Raises
        ------
        PropertyConversionException: If conversion is not possible
        """
        if prop_class is Resistivity:
            return Resistivity(prop.inverse_value, 'Automatically generated property from Conductivity.')
        if prop_class is DifferentialConductivity:
            if not prop.is_linear:
                logger.error("Conversion from \'Conductivity\' to \'DifferentialConductivity\' is only possible when "
                             "the conductivity is linear. The DifferentialConductivity has to be specified.")
                raise PropertyConversionException("Conversion from \'Conductivity\' to \'DifferentialConductivity\' is "
                                                  "only possible when the conductivity is linear")
            return DifferentialConductivity(prop.value, 'Automatically generated property from Conductivity.')
        raise PropertyConversionException(f"Conversion from \'Conductivity\' to {prop_class} not possible.")

    @staticmethod
    def from_density(prop: Density, prop_class):
        """Convert a MatProperty from a density

        Parameters
        ----------
        prop : Density
            Existing MatProperty.
        prop_class :
            Desired output property

        Returns
        -------
        new_prop :
            An object of the desired MatProperty or None, if not possible

        Raises
        ------
        PropertyConversionException: If conversion is not possible
        """
        raise PropertyConversionException(f"Conversion from \'Density\' to {prop_class} not possible.")

    @staticmethod
    def from_differential_conductivity(prop: DifferentialConductivity, prop_class):
        """Convert a MatProperty from a differential conductivity

        Parameters
        ----------
        prop : DifferentialConductivity
            Existing MatProperty.
        prop_class :
            Desired output property

        Returns
        -------
        new_prop :
            An object of the desired MatProperty or None, if not possible

        Raises
        ------
        PropertyConversionException: If conversion is not possible
        """
        raise PropertyConversionException(f"Conversion from \'DifferentialConductivity\' to {prop_class} not possible.")

    @staticmethod
    def from_differential_permittivity(prop: DifferentialPermittivity, prop_class):
        """Convert a MatProperty from a differential permittivity

        Parameters
        ----------
        prop : DifferentialPermittivity
            Existing MatProperty.
        prop_class :
            Desired output property

        Returns
        -------
        new_prop :
            An object of the desired MatProperty or None, if not possible

        Raises
        ------
        PropertyConversionException: If conversion is not possible
        """
        raise PropertyConversionException(f"Conversion from \'DifferentialPermittivity\' to {prop_class} not possible.")

    @staticmethod
    def from_permeability(prop: Permeability, prop_class):
        """Convert a MatProperty from a permeability

        Parameters
        ----------
        prop : Permeability
            Existing MatProperty.
        prop_class :
            Desired output property

        Returns
        -------
        new_prop :
            An object of the desired MatProperty or None, if not possible

        Raises
        ------
        PropertyConversionException: If conversion is not possible
        """
        if prop_class is Reluctivity:
            return Reluctivity(prop.inverse_value, 'Automatically generated property from Permeability.')
        raise PropertyConversionException(f"Conversion from \'Permeability\' to {prop_class} not possible.")

    @staticmethod
    def from_permittivity(prop: Permittivity, prop_class):
        """Convert a MatProperty from a permittivity

        Parameters
        ----------
        prop : Permittivity
            Existing MatProperty.
        prop_class :
            Desired output property

        Returns
        -------
        new_prop :
            An object of the desired MatProperty or None, if not possible

        Raises
        ------
        PropertyConversionException: If conversion is not possible
        """
        if prop_class is DifferentialPermittivity:
            if not prop.is_linear:
                logger.error("Conversion from \'Permittivity\' to \'DifferentialPermittivity\' is only possible when "
                             "the permittivity is linear. The DifferentialPermettivity has to be specified.")
                raise PropertyConversionException("Conversion from \'Permittivity\' to \'DifferentialPermittivity\' is "
                                                  "only possible when the permittivity is linear")
            return DifferentialPermittivity(prop.value, 'Automatically generated property from Permittivity.')
        raise PropertyConversionException(f"Conversion from \'Permittivity\' to {prop_class} not possible.")

    @staticmethod
    def from_reluctivity(prop: Reluctivity, prop_class):
        """Convert a MatProperty from a reluctivity

        Parameters
        ----------
        prop : Reluctivity
            Existing MatProperty.
        prop_class :
            Desired output property

        Returns
        -------
        new_prop :
            An object of the desired MatProperty or None, if not possible

        Raises
        ------
        PropertyConversionException: If conversion is not possible
        """
        if prop_class is Permeability:
            return Permeability(prop.inverse_value, 'Automatically generated property from Reluctivity.')
        raise PropertyConversionException(f"Conversion from \'Reluctivity\' to {prop_class} not possible.")

    @staticmethod
    def from_resistivity(prop: Resistivity, prop_class):
        """Convert a MatProperty from a resistivity

        Parameters
        ----------
        prop : Resistivity
            Existing MatProperty.
        prop_class :
            Desired output property

        Returns
        -------
        new_prop :
            An object of the desired MatProperty or None, if not possible

        Raises
        ------
        PropertyConversionException: If conversion is not possible
        """
        if prop_class is Conductivity:
            return Conductivity(prop.inverse_value, 'Automatically generated property from Resistivity.')
        raise PropertyConversionException(f"Conversion from \'Resistivity\' to {prop_class} not possible.")

    @staticmethod
    def from_thermal_conductivity(prop: ThermalConductivity, prop_class):
        """Convert a MatProperty from a thermal conductivity

        Parameters
        ----------
        prop : ThermalConductivity
            Existing MatProperty.
        prop_class :
            Desired output property

        Returns
        -------
        new_prop :
            An object of the desired MatProperty or None, if not possible

        Raises
        ------
        PropertyConversionException: If conversion is not possible
        """
        raise PropertyConversionException(f"Conversion from \'ThermalConductivity\' to {prop_class} not possible.")

    @staticmethod
    def from_volumetric_heat_capacity(prop: VolumetricHeatCapacity, prop_class):
        """Convert a MatProperty from a volumetric heat capacity

        Parameters
        ----------
        prop : VolumetricHeatCapacity
            Existing MatProperty.
        prop_class :
            Desired output property

        Returns
        -------
        new_prop :
            An object of the desired MatProperty or None, if not possible

        Raises
        ------
        PropertyConversionException: If conversion is not possible
        """
        raise PropertyConversionException(f"Conversion from \'VolumetricHeatCapacity\' to {prop_class} not possible.")


class Mat:
    """Class that represents a specific material."""

    __num_instances = 0  # Counts the number of instances of Mat
    __ID_count = 0  # Provides unique ID for the Mat

    __slots__ = ('_ID', 'name', '_mat_properties', 'property_converter')

    def __init__(self, name: str, *props: MatProperty):
        """
        Constructor of Mat

        Parameters
        ----------
        name : str
            Name of the Mat.
        *props : MatProperty
            A number of MatProperties, possibly empty.

        Returns
        -------
        None.

        """
        self._ID = Mat.__ID_count
        self.name = name
        self._mat_properties = []
        self.property_converter = PropertyConverter()

        if len(props) > 0:
            self.add_property(*props)

        # Increase both counters
        Mat.__num_instances += 1
        Mat.__ID_count += 1

    def __del__(self):
        """
        Delete the Mat

        Returns
        -------
        None.

        """
        Mat.__num_instances -= 1

    @property
    def ID(self):
        """
        Returns the ID of the Mat

        Returns
        -------
        int
            ID of the Mat.

        """
        return self._ID

    def is_linear(self, prop_class: Type[MatProperty] = None) -> bool:
        """Returns if this Mat is linear or not.

        Parameters
        ----------
        prop_class : Type[MatProperty], optional
            If given, it is only checked if the property of this type is linear. If not given, all properties are
            considered.

        Returns
        -------
        linear : bool
        """
        if prop_class is not None:
            return self.get_property(prop_class).is_linear
        return all((prop.is_linear for prop in self._mat_properties))

    def is_homogeneous(self, prop_class: Type[MatProperty] = None) -> bool:
        """Returns if this Mat is homogeneous or not.

        Parameters
        ----------
        prop_class : Type[MatProperty], optional
            If given, it is only checked if the property of this type is homogeneous. If not given, all properties are
            considered.

        Returns
        -------
        homogeneous : bool
        """
        if prop_class is not None:
            return self.get_property(prop_class).is_homogeneous
        return all((prop.is_homogeneous for prop in self._mat_properties))

    def is_isotropic(self, prop_class: Type[MatProperty] = None) -> bool:
        """Returns if this Mat is isotropic or not.

        Parameters
        ----------
        prop_class : Type[MatProperty], optional
            If given, it is only checked if the property of this type is isotropic. If not given, all properties are
            considered.

        Returns
        -------
        isotropic : bool
        """
        if prop_class is not None:
            return self.get_property(prop_class).is_isotrop
        return all((prop.is_isotrop for prop in self._mat_properties))

    def is_time_dependent(self, prop_class: Type[MatProperty] = None) -> bool:
        """Returns if this Mat is time dependent or not.

        Parameters
        ----------
        prop_class : Type[MatProperty], optional
            If given, it is only checked if the property of this type is time dependent. If not given, all properties
            are considered.

        Returns
        -------
        isotropic : bool
        """
        if prop_class is not None:
            return self.get_property(prop_class).is_time_dependent
        return any((prop.is_time_dependent for prop in self._mat_properties))

    def is_hysteretic(self, prop_class: Type[MatProperty] = None) -> bool:
        """Returns if this Mat is hysteretic or not.

        Parameters
        ----------
        prop_class : Type[MatProperty], optional
            If given, it is only checked if the property of this type is hysteretic. If not given, all properties are
            considered.

        Returns
        -------
        hysteretic : bool
        """
        if prop_class is not None:
            return self.get_property(prop_class).is_hysteretic
        return all((prop.is_hysteretic for prop in self._mat_properties))

    def update(self, name: str, val: Any, *, prop_class: Type[MatProperty] = None):
        """Update a value in the MatProperties.

        Parameters
        ----------
        name : str
            Name of the argument in the MatProperties to update.
        val : Any
            The new value that is written on the MatProperties
        prop_class : Type[MatProperty], optional
            If given, set the value only in the MatProperty of given type. If not given, set the value for all.
            Default is None
        """
        if prop_class is not None:
            mat_properties = [self.get_property(prop_class), ]
        else:
            mat_properties = self._mat_properties
        for mat_property in mat_properties:
            if mat_property is not None:
                mat_property.update(name, val)

    def add_property(self, *props: MatProperty):
        """
        Adds a number of MatProperties to the Mat

        Parameters
        ----------
        *props : MatProperty
            Properties of the Mat.

        Raises
        ------
        TypeError
            When an object in props is no instance of MatProperty.

        Returns
        -------
        None.

        """
        for prop in props:
            if not isinstance(prop, MatProperty):
                raise TypeError(
                    f"The property of type {str(type(prop))} is not a instance of class MatProperty ")
            self._mat_properties.append(prop)

    def has_original_property(self, prop_class) -> bool:
        """Returns True if a MatProperty of `prop_class` is directly defined on this Mat, False if not

        Parameters
        ----------
        prop_class :
            Class to check

        Returns
        -------
        has_property : bool
        """
        for prop in self._mat_properties:
            if isinstance(prop, prop_class):
                return True
        return False

    def _get_property_from_related(self, prop_class):
        """Create a new property of `prop_class` based on existing properties."""
        for prop in self._mat_properties:
            prop: MatProperty
            try:
                related_prop = self.property_converter.convert(prop, prop_class)
                return related_prop
            except PropertyConversionException:
                pass
        return None

    def get_property(self, prop_class) -> Union[Type[MatProperty], None]:
        """
        Returns the MatProperty of the class specified in prop_class

        Parameters
        ----------
        prop_class : MatProperty
            The class which MatProperty should be returned

        Returns
        -------
        MatProperty
            Object of class specified in prop_class.

        """
        for prop in self._mat_properties:
            if isinstance(prop, prop_class):
                return prop
        return self._get_property_from_related(prop_class)

    def delete_property(self, prop_class):
        """
        Deletes the MatProperty of prop_class from the Mat

        Parameters
        ----------
        prop_class :
            The class wich object should be deleted.

        Returns
        -------
        None.

        """
        for prop, k in zip(self._mat_properties, range(len(self._mat_properties))):
            if isinstance(prop, prop_class):
                del self._mat_properties[k]
                break
