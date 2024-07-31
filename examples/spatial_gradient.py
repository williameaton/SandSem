import matplotlib.pyplot as plt
from src.mesh.generate_mesh import generate_mesh
from src.solver.static_solver import StaticSolver
import numpy as np
from wetools.plotting import setup_we_mpl
setup_we_mpl()

# Create a basic mesh
m = generate_mesh(nproc_xi=16, nelem_rad=2, ngll=5, rmin=0.2, rmax=0.3)

# Create a static solver
s = StaticSolver(mesh=m)

# Generate a 1D variable (e.g. pressure)
s.generate_static_var(variable='pressure', ndim=1)

# Initialise an array for analytical solution
analytical = np.zeros((3,m.npoints))

# Compute pressure and analytical gradient of the pressure at each point
for i in range(m.npoints):
    x = m.gcoord_x[i]
    y = m.gcoord_y[i]
    z = m.gcoord_z[i]
    s.pressure[i] = (x**2)*z*y - 2*y*y + np.sin(z) + np.cos(30*x)

    analytical[0,i] = 2*x*y*z - 30*np.sin(30*x)
    analytical[1,i] = x*x*z - 4*y
    analytical[2,i] = x*x*y + np.cos(z)

# Compute the gradient using SEM
s.grad_pressure = s.compute_spatial_scalar_gradient(s.pressure)

# Compute raw difference
diff = s.grad_pressure - analytical


# Plot SEM, analytical and difference along surface (top boundary)
boundary = "top"

fig = plt.figure()
ax1 = fig.add_subplot(1, 3, 1, projection='3d')
ax2 = fig.add_subplot(1, 3, 2, projection='3d')
ax3 = fig.add_subplot(1, 3, 3, projection='3d')
axs = [ax1, ax2, ax3]

m.plot_variable_on_mesh_boundary(boundary_name=boundary, ax=ax1, variable=s.grad_pressure[0,:])
m.plot_variable_on_mesh_boundary(boundary_name=boundary, ax=ax2, variable=analytical[0,:])
m.plot_variable_on_mesh_boundary(boundary_name=boundary, ax=ax3, variable=diff[0,:] )


titles = ['SEM gradient\non surface', 'Analytical gradient\non surface', 'Difference']
for i in range(3):
    axs[i].set_title(titles[i])

plt.show()


