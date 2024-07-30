import matplotlib.pyplot as plt
from mesh.generate_mesh import generate_mesh
from solver.hyperbolic_solver import HyperbolicSolver

# Create a basic mesh
m = generate_mesh(nproc_xi=4, nelem_rad=2, ngll=3, rmin=0.2, rmax=0.3)
fig, ax = m.plot_mesh()

# Create a hyperbolic solver
s = HyperbolicSolver(mesh=m)

# Compute the mass matrix only once
#m.create_mass_matrix()

s.generate_variable_arrays('displacement')
s.zero_displacement()

m.compute_node_valency()

m.write_global_coords()

m.plot_global_variable(m.nvalency, fig,ax)

#s.plot_displacement()

#s.add_gaussian_source_in_element(amp, spatial_decay)


plt.show()
