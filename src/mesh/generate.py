from wetools.funcs import map_csb_to_xyz
from wetools.plotting import setup_we_mpl
import numpy as np
from mesh import Mesh
from element import Element
import matplotlib.pyplot as plt
from gll import mapping_deriv_3GLL_analytical
setup_we_mpl()

m = Mesh(nproc_xi=16, nelem_rad=3, chunk=1, ngll=5)


xi  = np.linspace(-np.pi/4, np.pi/4, m.nproc_xip1)
eta = np.linspace(-np.pi/4, np.pi/4, m.nproc_xip1)
r   = np.linspace(0.5, 1, m.nelem_radp1)

m.create_hex()

m.setup_ibool()

# Create element from its corners
for k in range(m.nelem_rad):
    for i in range(m.nproc_xi):
        for j in range(m.nproc_xi):
            elem = Element()
            elem.coord_x = np.zeros((m.ngll,m.ngll,m.ngll))
            elem.coord_y = np.zeros((m.ngll,m.ngll,m.ngll))
            elem.coord_z = np.zeros((m.ngll,m.ngll,m.ngll))

            xi1  = xi[i]
            xi2  = xi[i+1]
            eta1 = eta[j]
            eta2 = eta[j+1]
            r1   = r[k]
            r2   = r[k+1]

            # Corner nodes:
            for rz in [r1,r2]:
                elem.cnodes_xyz.append(map_csb_to_xyz(xi2, eta1, rz, m.chunk))
                elem.cnodes_xyz.append(map_csb_to_xyz(xi2, eta2, rz, m.chunk))
                elem.cnodes_xyz.append(map_csb_to_xyz(xi1, eta2, rz, m.chunk))
                elem.cnodes_xyz.append(map_csb_to_xyz(xi1, eta1, rz, m.chunk))

            # Interpolated nodes in CS basis:
            g = m.gll

            xi_int  = (xi2 + xi1)/2 + g*(xi2 - xi1)/2
            eta_int = (eta2 + eta1)/2 + g*(eta2 - eta1)/2
            r_int   = (r2 + r1)/2 + g*(r2 - r1)/2

            for rri in m.range_ngll:
                for eei in m.range_ngll:
                    for xxi in m.range_ngll:
                        x,y,z = map_csb_to_xyz(xi_int[xxi], eta_int[eei], r_int[rri], m.chunk)
                        elem.coord_x[xxi,eei,rri] = x
                        elem.coord_y[xxi,eei,rri] = y
                        elem.coord_z[xxi,eei,rri] = z
            elem.cnodes_xyz_to_np()
            m.add_element(elem)

print()

m.link_ibool_to_element()

m.set_homogenous_property('cijkl', 3e10)
m.set_homogenous_property('density', 2500)

m.setup_integration()

m.create_mass_matrix()






fig, ax = plt.subplots()
ax.plot(np.arange(m.npoints), m.M)
plt.show()


"""
fig = plt.figure()
fig.set_tight_layout(True)
ax = fig.add_subplot(projection='3d')

# Plot skeleton
for elem in m.elements:
    elem.plot_element(ax)

ax.set_xlim([-0.5,0.5])
ax.set_ylim([-1,1])
ax.set_zlim([-0.5,0.5])

plt.show()"""
