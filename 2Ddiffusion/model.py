import sys
import os.path
from input_utils import res_dir, K_names, terminated_flag, get_cmd_args_parser

# Read cmd line run params:
parser = get_cmd_args_parser()
args = parser.parse_args()

res_dir = res_dir.format(**vars(args))

layer_idx = args.layer_idx
npixels = args.npixels
active_vessels = args.active_vessels

if os.path.isfile(res_dir + K_names(terminated_flag)[1]):
    print('This run has already been executed: {}'.format(res_dir))
    sys.exit()

import numpy as np
from fenics import *
from fenics_np_utils import make_K_arr

PATH = '../VesselData/Meshes/'


print('Running on params:')
print(vars(args))
# Model parameters
dH_slice = 12e-6 # 12mum [Jakob]
C0      = args.C_0 # 15
Km      = args.K_m #1
K_vessel_surface \
        = args.K_0
D       = args.D * 1e-12 #* dH_slice #2000e-12
dt      = 20e-3   # time step
n_vessel_g = 20 # number of vessel groups used. This comes from the GMSH script make_gmesh_from_file..

if active_vessels%(100/n_vessel_g) != 0:
    raise ValueError('active vessels needs to be {} * x'.format(100/n_vessel_g))

mesh = Mesh()

# Full mesh:
with XDMFFile(PATH + "mesh_np={}_iz={}.xdmf".format(npixels, layer_idx)) as infile:
    infile.read(mesh)

# Mesh values/labels are stored here:
mvc = MeshValueCollection("size_t", mesh, 1)
with XDMFFile(PATH + "mf_np={}_iz={}.xdmf".format(npixels, layer_idx)) as infile:
    infile.read(mvc, "name_to_read")


mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)

# Create mesh and define function space
# TODO: What does the order of this basis mean?
V = FunctionSpace(mesh, 'CG', 1)


# Define boundary conditions.
# Vessels are labeld by numbers 10...n_vessel_g - 1

print('Active vessel groups')
vbc = []
for vessel_g in np.random.choice(np.arange(10, 10 + n_vessel_g), size=int(active_vessels/100*n_vessel_g), replace=False):
    print(vessel_g)
    vbc.append(DirichletBC(V, Constant(K_vessel_surface), mf, vessel_g))

# Define initial value
K_n = project(Constant(K_vessel_surface), V)

# Define variational problem
K = TrialFunction(V)
u = TestFunction(V)

F = u*K*dx + dt*dot(grad(u), grad(K))*D*dx \
        + u*C0*Km/(K_n+Km)**2*(K-K_n)*dt*dx - u*K_n*dx \
        + u*C0*K_n/(K_n + Km)*dt*dx
a, L = lhs(F), rhs(F)

# Create VTK file for saving solution
vtkfile = File(res_dir + 'K.pvd')

# Time-stepping
K = Function(V)
t = 0
i = 0
terminate = .0001


#for n in range(int(20/dt)):
while True:

    # Update current time using adative time stepping
    t += dt

    # Compute solution
    # TODO: make sure tha the von Neuman bcs are satified at the
    # domain boundary.
    solve(a == L, K, vbc)

    delta_K = (K_n.compute_vertex_values() - K.compute_vertex_values()) \
            /(K_n.compute_vertex_values() + Km)
    print('{:.05f}'.format(delta_K.max()),
          delta_K.argmax(),
          (delta_K > .01).sum(),
          len(delta_K), dt)


    # Update previous solution
    K_n.assign(K)

    # Start saving results and plotting:
    if i%500==0:
        print('Time: {:.03f}ms'.format(t*1000))
        #vtkfile << (K, t)
        #data = make_K_arr(mesh, K_n, t, res_dir)

    if delta_K.max() < terminate:
        make_K_arr(mesh, K_n, terminated_flag, res_dir)
        print('Terminated at t={:.02f}s'.format(t))
        break

    i += 1
