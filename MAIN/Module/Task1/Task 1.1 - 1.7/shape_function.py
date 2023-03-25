import numpy as np
from dataclasses import dataclass

# Definition of shape function N = a+bx+cy
'''
einfache data Klasse für Shapefunction
'''
@dataclass
class ShapeFunction_N:

    def __init__(self, depth):
        self.depth = depth

    def calc_coeff(self, msh):
        x1 = msh.node[msh.elem_to_node, 0]  # x coordinates of the nodes for each element (num_nodes x 1)
        y1 = msh.node[msh.elem_to_node, 1]  # y coordinates of the nodes for each element (num_nodes x 1)

        x2 = np.roll(x1, -1, axis=1)
        y2 = np.roll(y1, -1, axis=1)
        x3 = np.roll(x1, -2, axis=1)
        y3 = np.roll(y1, -2, axis=1)

        # für jeden Punkt die Nodal(edge) function
        self.a = x2 * y3 - x3 * y2
        self.b = y2 - y3
        self.c = x3 - x2
        self.element_area = np.mean(((x2 * y3 - x3 * y2) + (y2 - y3) * x1 + (x3 - x2) * y1) / 2, 1)

    def get_coeff(self):
        return self.a, self.b, self.c, self.element_area




