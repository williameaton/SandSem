from wetools.funcs import map_csb_to_xyz
import numpy as np
from classes import Mesh, Node, Element
import matplotlib.pyplot as plt


m = Mesh()

# Specify boundaries
ngll = 5
m.chunk = 1
m.maxdepth = 1000 # km

nproc_xi   = 6
nproc_xip1 = nproc_xi+1

xi  = np.linspace(-np.pi/4, np.pi/4, nproc_xip1)
eta = np.linspace(-np.pi/4, np.pi/4, nproc_xip1)
r   = np.linspace(0.5, 1, nproc_xip1)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')


for i in range(nproc_xi):
    for j in range(nproc_xi):
        # Create element from its corners

        elem = Element()
        for rz in range(2):

            n = Node()
            n.x, n.y, n.z = map_csb_to_xyz(xi[i], eta[j], r[rz], m.chunk)
            n.ngll = ngll
            n.id   = 1 + ngll*(ngll-1)*rz
            elem.add_corner_node(n)


            n = Node()
            n.x, n.y, n.z = map_csb_to_xyz(xi[i+1], eta[j], r[rz], m.chunk)
            n.ngll = ngll
            n.id   = ngll  + ngll*(ngll-1)*rz
            elem.add_corner_node(n)

            n = Node()
            n.x, n.y, n.z = map_csb_to_xyz(xi[i+1], eta[j+1], r[rz], m.chunk)
            n.ngll = ngll
            n.id   = ngll**3 - (ngll-1)*ngll  + ngll*(ngll-1)*rz

            elem.add_corner_node(n)

            n = Node()
            n.x, n.y, n.z = map_csb_to_xyz(xi[i+1], eta[j+1], r[rz], m.chunk)
            n.ngll = ngll
            n.id   = ngll**3 - ngll**2  + 1  + ngll*(ngll-1)*rz
            elem.add_corner_node(n)

        m.add_element(elem)

ax.set_xlim([-1,1])
ax.set_ylim([-1,1])
ax.set_zlim([-1,1])

plt.show()
