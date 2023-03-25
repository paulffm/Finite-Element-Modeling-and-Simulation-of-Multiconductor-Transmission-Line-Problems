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
        x1 = msh.node[msh.elem_to_node, 0]  # x coordinates of the nodes for each element (num_nodes x 3)
        y1 = msh.node[msh.elem_to_node, 1]  # y coordinates of the nodes for each element (num_nodes x 3)


        a = np.zeros((x1.shape[0], x1.shape[1]))
        b = np.zeros((x1.shape[0], x1.shape[1]))
        c = np.zeros((x1.shape[0], x1.shape[1]))

        "ai = xj * yk - xk * yj"
        "bi = yj - yk"
        "ci = xk - xj"


        for i in range((x1.shape[1])):
            if i == 0:
                'i=0, j=1, k=2'
                a[:, i] = x1[:, 1] - y1[:, 2] - x1[:, 2] - y1[:, 1]
                b[:, i] = y1[:, 1] - y1[:, 2]
                c[:, i] = x1[:, 2] - x1[:, 1]

            if i == 1:
                'i=1, j=2, k=1'
                a[:, i] = x1[:, 2] - y1[:, 0] - x1[:, 0] - y1[:, 2]
                b[:, i] = y1[:, 2] - y1[:, 0]
                c[:, i] = x1[:, 0] - x1[:, 2]

            if i == 2:
                'i=2, j=0, k=1'
                a[:, i] = x1[:, 0] - y1[:, 1] - x1[:, 1] - y1[:, 0]
                b[:, i] = y1[:, 0] - y1[:, 1]
                c[:, i] = x1[:, 1] - x1[:, 0]

        self.a = a
        self.b = b
        self.c = c

        "calculates the element area with the cross product 1/2*|a_vec x b_vec|"
        ax = x1[:, 1] - x1[:, 0]
        ay = y1[:, 1] - y1[:, 0]
        bx = x1[:, 2] - x1[:, 0]
        by = y1[:, 2] - y1[:, 0]
        self.element_area = 0.5 * np.abs(ax * by - ay * bx)


    def get_coeff(self):
        return self.a, self.b, self.c, self.element_area




