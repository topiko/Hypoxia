from fenics import *
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from input_utils import K_names

def make_K_arr(mesh, K_n, t, res_dir):
    # This is to verify sanity using pyplot

    # TODO: Better way to count cells!???
    ncells = 0
    for cell in cells(mesh):
        ncells += 1
    data = np.zeros((ncells*3, 6))

    i = 0
    for cell in cells(mesh):
        coord_arr = np.array(cell.get_vertex_coordinates()).reshape(3,3)
        x = coord_arr[:, 0]
        y = coord_arr[:, 1]
        z = coord_arr[:, 2]

        for coord in coord_arr:
            data[i,0] = cell.index()
            data[i,1:4] = coord
            data[i,4] = K_n(Point(*coord))
            data[i,5] = cell.volume()
            i += 1

    fname = res_dir + K_names(t)[1]

    np.save(fname, data)

    print('saving data in: "{}"'.format(fname))

    return data

