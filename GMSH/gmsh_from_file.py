# ------------------------------------------------------------------------------
#
#  Gmsh Python tutorial 1
#
#  Elementary entities (points, curves, surfaces), physical groups
#  (points, curves, surfaces)
#
# ------------------------------------------------------------------------------

# The Python API is entirely defined in the `gmsh.py' module (which contains the
# full documentation of all the functions in the API):
import gmsh
import sys
import numpy as np

# Before using any functions in the Python API, Gmsh must be initialized:
gmsh.initialize()

# By default Gmsh will not print out any messages: in order to output messages
# on the terminal, just set the "General.Terminal" option to 1:
gmsh.option.setNumber("General.Terminal", 1)

# Next we add a new model named "t1" (if gmsh.model.add() is not called a new
# unnamed model will be created on the fly, if necessary):
gmsh.model.add("wo_vessels")

# The Python API provides direct access to each supported geometry kernel. The
# built-in kernel is used in this first tutorial: the corresponding API
# functions have the `gmsh.model.geo' prefix.

# The first type of `elementary entity' in Gmsh is a `Point'. To create a point
# with the built-in geometry kernel, the Python API function is
# gmsh.model.geo.addPoint():
# - the first 3 arguments are the point coordinates (x, y, z)
# - the next (optional) argument is the target mesh size (the "characteristic
#   length") close to the point
# - the last (optional) argument is the point tag (a stricly positive integer
#   that uniquely identifies the point)
lc = 2e-5
lc2 = 2e-5
#xmin = 0
#xmax = 1
#ymin = 0
#ymax = 1
#xmin_h = 0.2
#xmax_h = 0.3
#ymin_h = 0.25
#ymax_h = 0.3

PATH = '../VesselData/'

if len(sys.argv) > 1:
    layer_idx = sys.argv[1]
    npixels = sys.argv[2]
    if len(sys.argv) == 4:
        show = sys.argv[3]
    else:
        show = False
else:
    print('Using default params: layer_idx = 1, npixels = 100')
    layer_idx = 1
    npixels = 100

print('Generating for layer {}'.format(layer_idx))
vessel_dict_init = np.load(PATH + 'vessels_dict_np={}_iz={}.npy'.format(npixels,
                                                                        layer_idx),
                           allow_pickle=True).item()

add_space = 5e-6 #.00001
xmin = -add_space
ymin = -add_space
xmax = vessel_dict_init['W'] + add_space
ymax = vessel_dict_init['H'] + add_space

print(xmin, xmax, ymin, ymax)
min_distance = 2.5e-6

# Filtre all but vessels off..
vessel_dict = {}
count = 0
for key, val in vessel_dict_init.items():
    if key.startswith('vessel'):
        vessel_dict[key] = val
        vessel_dict[key]['point_tags'] = []
        vessel_dict[key]['line_tags'] = []
        count += 1

# The distribution of the mesh element sizes will be obtained by interpolation
# of these characteristic lengths throughout the geometry. Another method to
# specify characteristic lengths is to use general mesh size Fields (see
# `t10.py'). A particular case is the use of a background mesh (see `t7.py').
#
# If no target mesh size of provided, a default uniform coarse size will be used
# for the model, based on the overall model size.
#
# We can then define some additional points. All points should have different
# tags:
gmsh.model.geo.addPoint(xmin, ymin, 0, lc, 1)
gmsh.model.geo.addPoint(xmax, ymin,  0, lc, 2)
gmsh.model.geo.addPoint(xmax, ymax, 0, lc, 3)
gmsh.model.geo.addPoint(xmin, ymax, 0, lc, 4)

for vessel in vessel_dict:
    for x, y in vessel_dict[vessel]['edge']:
        if (x in [xmin, xmax]) or (y in [ymin, ymax]):
            raise ValueError('Edge violated')
        p = gmsh.model.geo.addPoint(x, y, 0, lc2)
        vessel_dict[vessel]['point_tags'].append(p)

# Curves are Gmsh's second type of elementery entities, and, amongst curves,
# straight lines are the simplest. The API to create straight line segments with
# the built-in kernel follows the same conventions: the first 2 arguments are
# point tags (the start and end points of the line), and the last (optional one)
# is the line tag.
#
# In the commands below, for example, the line 1 starts at point 1 and ends at
# point 2.
#
# Note that curve tags are separate from point tags - hence we can reuse tag `1'
# for our first curve. And as a general rule, elementary entity tags in Gmsh
# have to be unique per geometrical dimension.
gmsh.model.geo.addLine(1, 2, 1)
gmsh.model.geo.addLine(2, 3, 2)
gmsh.model.geo.addLine(3, 4, 3)
gmsh.model.geo.addLine(4, 1, 4)

#vessel_line_tags = []

for vessel in vessel_dict:
    for tag in vessel_dict[vessel]['point_tags'][:-1]: #vessel_point_tags[:-1]:
        p = gmsh.model.geo.addLine(tag, tag+1)
        vessel_dict[vessel]['line_tags'].append(p)
        #vessel_line_tags.append(p)
    p = gmsh.model.geo.addLine(vessel_dict[vessel]['point_tags'][-1],
                               vessel_dict[vessel]['point_tags'][0])
    vessel_dict[vessel]['line_tags'].append(p)

#gmsh.model.geo.addLine(6, 7, 6)
#gmsh.model.geo.addLine(7, 8, 7)
#gmsh.model.geo.addLine(8, 5, 8)


# The third elementary entity is the surface. In order to define a simple
# rectangular surface from the four curves defined above, a curve loop has first
# to be defined. A curve loop is defined by an ordered list of connected curves,
# a sign being associated with each curve (depending on the orientation of the
# curve to form a loop). The API function to create curve loops takes a list
# of integers as first argument, and the curve loop tag (which must ne unique
# amongst curve loops) as the second (optional) argument:

vessel_tags = []
for vessel in vessel_dict:
    vessel_line_tags = vessel_dict[vessel]['line_tags']
    p = gmsh.model.geo.addCurveLoop(vessel_line_tags)
    vessel_dict[vessel]['loop_tag'] = p
    vessel_dict[vessel]['surface_tag'] = p
    gmsh.model.geo.addPlaneSurface([p], p)
    vessel_tags.append(p)

gmsh.model.geo.addCurveLoop([4, 1, 2, 3], p+1)

# We can then define the surface as a list of curve loops (only one here,
# representing the external contour, since there are no holes--see `t4.py' for
# an example of a surface with a hole):

#gmsh.model.geo.addPlaneSurface([1], 1)
gmsh.model.geo.addPlaneSurface([p+1, *vessel_tags], p+1)

# At this level, Gmsh knows everything to display the rectangular surface 1 and
# to mesh it. An optional step is needed if we want to group elementary
# geometrical entities into more meaningful groups, e.g. to define some
# mathematical ("domain", "boundary"), functional ("left wing", "fuselage") or
# material ("steel", "carbon") properties.
#
# Such groups are called "Physical Groups" in Gmsh. By default, if physical
# groups are defined, Gmsh will export in output files only mesh elements that
# belong to at least one physical group. (To force Gmsh to save all elements,
# whether they belong to physical groups or not, set the `Mesh.SaveAll' option
# to 1.) Physical groups are also identified by tags, i.e. stricly positive
# integers, that should be unique per dimension (0D, 1D, 2D or 3D). Physical
# groups can also be given names.
#
# Here we define a physical curve that groups the left, bottom and right curves
# in a single group (with prescribed tag 5); and a physical surface with name
# "My surface" (with an automatic tag) containing the geometrical surface 1:

# Build vessel groups:
# The tags are later used in fenics to define bcs!
n_vessel_groups = 20
success = False
while not success:
    vessel_line_tags = {i:[] for i in range(n_vessel_groups)}
    for vessel in vessel_dict:
        # Each vessel is randomly allocated to a group:
        i = np.random.randint(0, n_vessel_groups)
        vessel_line_tags[i] += vessel_dict[vessel]['line_tags']

    # Check that ech group has at least single vessel if not redo...
    success = True
    for i in range(n_vessel_groups):
        if len(vessel_line_tags[i]) == 0:
            break
            success = False

    if not success:
        continue

    print('passed')
    for i in range(n_vessel_groups):
        # Build groups ad tags:
        vessel_group_idx = 10 + i
        gmsh.model.addPhysicalGroup(1, vessel_line_tags[i], vessel_group_idx)
        gmsh.model.setPhysicalName(1, vessel_group_idx, "vessels_{:02d}".format(vessel_group_idx))

gmsh.model.addPhysicalGroup(1, [1, 2, 3, 4], 5)
gmsh.model.setPhysicalName(1, 5, "exterior")

gmsh.model.addPhysicalGroup(2, [p+1], 7)
gmsh.model.setPhysicalName(1, 7, "tissue")


# Before it can be meshed, the internal CAD representation must be synchronized
# with the Gmsh model, which will create the relevant Gmsh data structures. This
# is achieved by the gmsh.model.geo.synchronize() API call for the built-in
# geometry kernel. Synchronizations can be called at any time, but they involve
# a non trivial amount of processing; so while you could synchronize the
# internal CAD data after every CAD command, it is usually better to minimize
# the number of synchronization points.
gmsh.model.geo.synchronize()

gmsh.model.setColor([(1, 1), (1, 2), (1, 3), (1, 4)], 0, 0, 255)
gmsh.model.setColor([(2, p['surface_tag']) for _, p in vessel_dict.items()], 255, 0, 0)

# We can then generate a 2D mesh...
gmsh.model.mesh.generate(2)

# Remember that by default, if physical groups are defined, Gmsh will export in
# the output mesh file only those elements that belong to at least one physical
# group. To force Gmsh to save all elements, you can use
#
#gmsh.option.setNumber("Mesh.SaveAll", 1)

# ... and save it to disk
gmsh.write("../VesselData/Meshes/w_vessels_np={}_iz={}.msh".format(npixels, layer_idx))

if show:
    gmsh.fltk.run()

# This should be called when you are done using the Gmsh Python API:
gmsh.finalize()
