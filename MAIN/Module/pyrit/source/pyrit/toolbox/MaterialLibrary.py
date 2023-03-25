# coding=utf-8
"""Implementation of a material library

.. sectionauthor:: Bundschuh, Ruppert
"""
import numpy as np
from scipy.constants import mu_0, epsilon_0
from pyrit.material.Permittivity import Permittivity
from pyrit.material.Permeability import Permeability
from pyrit.material.Conductivity import Conductivity
from pyrit.material.Density import Density
from pyrit.material.VolumetricHeatCapacity import VolumetricHeatCapacity
from pyrit.material.ThermalConductivity import ThermalConductivity
from pyrit.material.Reluctivity import Reluctivity
from pyrit.material.Mat import Mat

#: Tuple[str]: Tuple with allowed identifiers for :py:func:`get_material`
present_materials = ('air', 'copper', 'iron')


def get_material(identifier: str) -> Mat:
    """Return an instance of :py:class:`~pyrit.material.Mat.Mat` with the properties of the specified material.

    Parameters
    ----------
    identifier : str
        Identifier of the material. See the variable :py:const:`~pyrit.toolbox.MaterialLibrary.present_materials`

    Returns
    -------
    Mat : Mat
        An instance of :py:class:`~pyrit.material.Mat.Mat` with all properties set.
    """
    if identifier not in present_materials:
        raise KeyError(f"There is no material for given identifier {identifier}.")

    properties = {}

    if identifier == 'air':
        properties = {'conductivity': 0, 'density': 1.189, 'volumetricheatcapacity': 1.006,
                      'thermalconductivity': 25.87,
                      'permittivity': 1.00059, 'permeability': 1.00000037}

    if identifier == 'copper':
        properties = {'conductivity': 58e6, 'density': 8920, 'volumetricheatcapacity': 385,
                      'thermalconductivity': 400,
                      'permeability': 0.999994}

    if identifier == 'iron':
        properties = {'conductivity': 11.7e6, 'density': 7847, 'volumetricheatcapacity': 449,
                      'thermalconductivity': 80,
                      'permeability': 5000}

    values = {'conductivity': np.NAN, 'density': np.NAN, 'volumetricheatcapacity': np.NAN,
              'thermalconductivity': np.NAN,
              'permittivity': np.NAN, 'permeability': np.NAN}
    for key, val in properties.items():
        values[key] = val

    return Mat(identifier, Conductivity(values['conductivity']), Density(values['density']),
               VolumetricHeatCapacity(values['volumetricheatcapacity']),
               ThermalConductivity(values['thermalconductivity']),
               Permittivity(epsilon_0 * values['permittivity']), Permeability(mu_0 * values['permeability']),
               Reluctivity(1 / (mu_0 * values['permeability'])))
