# coding=utf-8
"""Combines different json files with integration schemes to a single file.

The single json files were generated with the matlab script 'write_json.m'.
The schemes can be downloaded from here: https://www.math.unipd.it/~alvise/SETS_CUBATURE_TRIANGLE/rules_triangle.html

.. sectionauthor:: Bundschuh
"""
# pylint: disable=unspecified-encoding

import json
import numpy as np


def main():
    """Combine different json files with integration schemes to a single file."""
    d = {}
    for k in range(1, 21):
        with open(f'dunavant_{k}.json', 'r') as f:
            data = json.load(f)

        xyz = np.array(data['xyz'])
        if xyz.ndim == 1:
            points = xyz[0:2]
            weights = xyz[2]
            d[f'{k}'] = {'points': points.tolist(), 'weights': weights}
        else:
            points = xyz[:, 0:2]
            weights = xyz[:, 2]

            d[f'{k}'] = {'points': points.tolist(), 'weights': weights.tolist()}

    with open('dunavant_all.json', 'w') as json_file:
        json.dump(d, json_file, indent=4, sort_keys=True)


if __name__ == '__main__':
    main()
