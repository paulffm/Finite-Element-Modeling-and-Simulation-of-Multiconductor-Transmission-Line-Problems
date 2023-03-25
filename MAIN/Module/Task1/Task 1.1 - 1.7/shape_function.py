import numpy as np
from dataclasses import dataclass

# Definition of shape function N = a+bx+cy
'''
einfache data Klasse für Shapefunction
'''
@dataclass
class ShapeFunction_N:
    depth: float
    element_area: np.ndarray
    a: np.ndarray
    b: np.ndarray
    c: np.ndarray


