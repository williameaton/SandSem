from wetools.funcs import map_csb_to_xyz
from wetools.plotting import setup_we_mpl
import numpy as np
from mesh import Mesh
from element import Element
import matplotlib.pyplot as plt
setup_we_mpl()

m = Mesh(nproc_xi=8, chunk=1, ngll=5)
m.maxdepth = 1000 # km

xi  = np.linspace(-np.pi/4, np.pi/4, m.nproc_xip1)
eta = np.linspace(-np.pi/4, np.pi/4, m.nproc_xip1)
r   = np.linspace(0.5, 1, m.nproc_xip1)


m.create_hex()

# Create element from its corners
for i in range(m.nproc_xi):
    for j in range(m.nproc_xi):
        elem = Element()
        for rz in range(2):
            elem.cnodes_xyz.append(map_csb_to_xyz(xi[i+1], eta[j],   r[rz], m.chunk))
            elem.cnodes_xyz.append(map_csb_to_xyz(xi[i+1], eta[j+1], r[rz], m.chunk))
            elem.cnodes_xyz.append(map_csb_to_xyz(xi[i],   eta[j+1], r[rz], m.chunk))
            elem.cnodes_xyz.append(map_csb_to_xyz(xi[i],   eta[j],   r[rz], m.chunk))

        elem.cnodes_xyz_to_np()
        m.add_element(elem)



# Add GLL points to elements
m.setup_gll()





fig = plt.figure()
fig.set_tight_layout(True)
ax = fig.add_subplot(projection='3d')


# Plot skeleton
for elem in m.elements:
    elem.plot_element(ax)


ax.set_xlim([-0.5,0.5])
ax.set_ylim([-1,1])
ax.set_zlim([-0.5,0.5])

plt.show()
