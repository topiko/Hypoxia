import meshio
import numpy as np
import sys
#from dolfin import *

iz = int(sys.argv[1])
npixels = int(sys.argv[2])

PATH = '../VesselData/Meshes/'
msh = meshio.read(PATH + "w_vessels_np={}_iz={}.msh".format(npixels, iz))

line_cells = []
for cell in msh.cells:
    if cell.type == "triangle":
        triangle_cells = cell.data
    elif  cell.type == "line":
        if len(line_cells) == 0:
            line_cells = cell.data
        else:
            line_cells = np.vstack([line_cells, cell.data])

line_data = []
# gmesh stores the various label to physical:
for key in msh.cell_data_dict["gmsh:physical"].keys():
    if key == "line":
        if len(line_data) == 0:
            line_data = msh.cell_data_dict["gmsh:physical"][key]
        else:
            line_data = np.vstack([line_data, msh.cell_data_dict["gmsh:physical"][key]])
    elif key == "triangle":
        triangle_data = msh.cell_data_dict["gmsh:physical"][key]

triangle_mesh = meshio.Mesh(points=msh.points,
                            cells={"triangle": triangle_cells})
line_mesh = meshio.Mesh(points=msh.points,
                        cells=[("line", line_cells)],
                        cell_data={"name_to_read":[line_data]})

meshio.write(PATH + "mesh_np={}_iz={}.xdmf".format(npixels, iz), triangle_mesh)
meshio.xdmf.write(PATH + "mf_np={}_iz={}.xdmf".format(npixels, iz), line_mesh)


