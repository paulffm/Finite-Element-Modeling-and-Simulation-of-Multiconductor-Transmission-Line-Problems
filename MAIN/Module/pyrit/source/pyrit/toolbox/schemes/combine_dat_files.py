# coding=utf-8
"""Combining dat files with cubature schemes to a json file.

The files can be downloaded from here: https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.6313

.. sectionauthor:: Bundschuh
"""
# pylint: disable=unspecified-encoding

import json
import numpy as np

file_names = ["NME_6313_cubature_tetra_p19_n390.dat",
              "NME_6313_cubature_tetra_p19_n392.dat",
              "NME_6313_cubature_tetra_p20_n448.dat",
              "NME_6313_cubature_tetra_p2_n4.dat",
              "NME_6313_cubature_tetra_p3_n6.dat",
              "NME_6313_cubature_tetra_p4_n11.dat",
              "NME_6313_cubature_tetra_p5_n14.dat",
              "NME_6313_cubature_tetra_p6_n23.dat",
              "NME_6313_cubature_tetra_p7_n31.dat",
              "NME_6313_cubature_tetra_p8_n44.dat",
              "NME_6313_cubature_tetra_p9_n57.dat",
              "NME_6313_cubature_tetra_p10_n74.dat",
              "NME_6313_cubature_tetra_p11_n94.dat",
              "NME_6313_cubature_tetra_p12_n117.dat",
              "NME_6313_cubature_tetra_p13_n144.dat",
              "NME_6313_cubature_tetra_p14_n175.dat",
              "NME_6313_cubature_tetra_p15_n207.dat",
              "NME_6313_cubature_tetra_p16_n247.dat",
              "NME_6313_cubature_tetra_p17_n288.dat",
              "NME_6313_cubature_tetra_p18_n338.dat"]


def main():
    """Combining dat files with cubature schemes to a json file"""
    data = {'1': {'points': [0.25, 0.25, 0.25], 'weights': [0.16666666666666666, ]}}
    for file_name in file_names:
        with open(file_name, 'r') as f:
            lines = f.readlines()
            points_and_weights = np.fromstring(''.join(lines), sep=' ').reshape((-1, 4))
            degree = file_name.split('_')[4][1:]
            data[degree] = {'points': points_and_weights[:, :3].tolist(),
                            'weights': (1 / 6 * points_and_weights[:, -1]).tolist()}

    with open('jaskowiec_sukumar_all.json', 'w') as json_file:
        json.dump(data, json_file, indent=4, sort_keys=True)


if __name__ == '__main__':
    main()
