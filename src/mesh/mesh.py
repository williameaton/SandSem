import matplotlib.pyplot as plt
import numpy as np
from wetools.plotting import setup_we_mpl
from wetools.funcs import map_csb_to_xyz
from src.mesh.gll import gll, lagrange1st
from src.mesh.hex import Hex
from src.mesh.element import Element


class Mesh():

    def __init__(self, nproc_xi, nelem_rad, chunk=1, ngll=5):

        self.nelem_rad   = nelem_rad
        self.nelem_radp1 = nelem_rad + 1
        self.nproc_xi   = nproc_xi
        self.nproc_xip1 = nproc_xi + 1
        self.ngll       = ngll
        self.chunk      = chunk

        self.lats = None
        self.lons = None
        self.elements   = []
        self.nelmts = 0

        self.homogeneous = {}

    def add_element(self,e):
        # Add element and link to mesh
        e.hex  = self.hex
        e.ngll = self.ngll
        e.mesh = self

        self.elements.append(e)
        self.nelmts += 1

    def setup_ibool(self):
        # get GLL points and weights:
        n = self.ngll
        self.gll, self.wgll = gll(n-1)
        self.ldash = lagrange1st(n - 1)
        self.range_ngll = range(n)

        # Total number of independent GLL points
        self.npoints_xi  = self.nproc_xip1 + self.nproc_xi * (n - 2)
        self.nrad_layers = self.nelem_rad*(n-1) + 1
        self.npoints     = (self.npoints_xi**2) * self.nrad_layers
        self.nelmts_hyp  = self.nproc_xi**2 * self.nelem_rad
        self.ibool       = np.zeros((n, n, n, self.nelmts_hyp), int)

        iglob = np.arange(self.npoints)
        self.iglob = iglob.reshape((self.npoints_xi, self.npoints_xi, self.nrad_layers), order='F')  # only 1 element deep

        ielmt = 0
        for k in range(self.nelem_rad):
            for i in range(self.nproc_xi):
                for j in range(self.nproc_xi):
                    for izeta in range(n):
                        for ieta in range(n):
                            for ixi in range(n):

                                iind = i*(n-1) + ixi
                                jind = j*(n-1) + ieta
                                kind = k*(n-1) + izeta

                                aa = self.iglob[iind,  jind,  kind]
                                self.ibool[ixi, ieta, izeta, ielmt] = aa
                    ielmt += 1


    def create_hex(self):
        self.hex = Hex(self.ngll)


    def link_ibool_to_element(self):

        for iel in range(self.nelmts):
            self.elements[iel].ibool = self.ibool[:,:,:,iel]


    def set_homogenous_property(self, type, value):
        # Sets a property of the mesh to be homogenous in all elements:
        if type == 'density':
            arr = np.zeros((self.ngll, self.ngll, self.ngll)) + value
        elif type == 'cijkl':
            arr = np.zeros((3,3,3,3, self.ngll, self.ngll, self.ngll)) + value
        else:
            raise ValueError('type not implemented yet')

        # Point elements to this array:
        for el in self.elements:
            el.prop[type] = arr


    def setup_elements(self):
        xi  = np.linspace(-np.pi / 4, np.pi / 4, self.nproc_xip1)
        eta = np.linspace(-np.pi / 4, np.pi / 4, self.nproc_xip1)
        r   = np.linspace(self.rmin, self.rmax, self.nelem_radp1)

        n = self.ngll
        # Create element from its corners
        for k in range(self.nelem_rad):
            for i in range(self.nproc_xi):
                for j in range(self.nproc_xi):
                    elem = Element()
                    elem.coord_x = np.zeros((n, n, n))
                    elem.coord_y = np.zeros((n, n, n))
                    elem.coord_z = np.zeros((n, n, n))

                    xi1 = xi[i]
                    xi2 = xi[i + 1]
                    eta1 = eta[j]
                    eta2 = eta[j + 1]
                    r1 = r[k]
                    r2 = r[k + 1]

                    # Corner nodes:
                    for rz in [r1, r2]:
                        elem.cnodes_xyz.append(map_csb_to_xyz(xi2, eta1, rz, self.chunk))
                        elem.cnodes_xyz.append(map_csb_to_xyz(xi2, eta2, rz, self.chunk))
                        elem.cnodes_xyz.append(map_csb_to_xyz(xi1, eta2, rz, self.chunk))
                        elem.cnodes_xyz.append(map_csb_to_xyz(xi1, eta1, rz, self.chunk))

                    # Interpolated nodes in CS basis:
                    g = self.gll

                    xi_int = (xi2 + xi1) / 2 + g * (xi2 - xi1) / 2
                    eta_int = (eta2 + eta1) / 2 + g * (eta2 - eta1) / 2
                    r_int = (r2 + r1) / 2 + g * (r2 - r1) / 2

                    for rri in self.range_ngll:
                        for eei in self.range_ngll:
                            for xxi in self.range_ngll:
                                x, y, z = map_csb_to_xyz(xi_int[xxi], eta_int[eei], r_int[rri], self.chunk)
                                elem.coord_x[xxi, eei, rri] = x
                                elem.coord_y[xxi, eei, rri] = y
                                elem.coord_z[xxi, eei, rri] = z
                    elem.cnodes_xyz_to_np()
                    self.add_element(elem)


    def write_global_coords(self):
        # Takes the element-wise defined coordinates and creates global arrays:
        self.gcoord_x = np.zeros((self.npoints))
        self.gcoord_y = np.zeros((self.npoints))
        self.gcoord_z = np.zeros((self.npoints))

        for ielmt in self.elements:
            ib = ielmt.ibool
            for a in self.range_ngll:
                for b in self.range_ngll:
                    for g in self.range_ngll:

                        iib = ib[a,b,g]
                        self.gcoord_x[iib] = ielmt.coord_x[a,b,g]
                        self.gcoord_y[iib] = ielmt.coord_y[a,b,g]
                        self.gcoord_z[iib] = ielmt.coord_z[a,b,g]


    def setup_integration(self):
        self.compute_jacobians()
        self.compute_weights()


    def compute_jacobians(self):
        for e in self.elements:
            e.compute_jacobian()


    def compute_weights(self):
        self.w3D = np.zeros((self.ngll, self.ngll, self.ngll))
        for a in self.range_ngll:
            for b in self.range_ngll:
                for g in self.range_ngll:
                    self.w3D[a,b,g] = self.wgll[a]*self.wgll[b]*self.wgll[g]
        # Point to it in elements
        # Main reason here is we could override it in the future if required
        for el in self.elements:
            el.w3D = self.w3D


    def create_mass_matrix(self):
        print('generating mass matrix...', end="")

        assert(self.nelmts == self.nelmts_hyp)
        self.M = np.zeros((self.npoints )) # Diagonal mass matrix:

        for i_elmt in range(self.nelmts):
            e = self.elements[i_elmt]
            ib = e.ibool
            rho = e.prop['density']

            for a in self.range_ngll:
                for b in self.range_ngll:
                    for g in self.range_ngll:

                        self.M[ib[a,b,g]] += rho[a,b,g] * e.jac3D[a,b,g] * e.w3D[a,b,g]

        print('done!')



    def compute_node_valency(self):
        self.nvalency = np.zeros((self.npoints))
        self.map_local_to_global(local=np.ones((self.ngll,
                                                self.ngll,
                                                self.ngll,
                                                self.nelmts)),
                                 glob=self.nvalency,
                                 avg_elmts=False)


    def map_local_to_global(self, local, glob=None, avg_elmts=False):
        # maps a local (elemental) array to a global array
        # values are averaged across the elemental values using the valency
        if (glob==None).any():
            glob = np.zeros((self.npoints))
        else:
            glob[:] = 0

        for ielmt in range(self.nelmts):
            ib = self.elements[ielmt].ibool
            for a in self.range_ngll:
                for b in self.range_ngll:
                    for g in self.range_ngll:
                        glob[ib[a,b,g]] += local[a,b,g,ielmt]

        # Determine whether to average
        if avg_elmts:
            return glob/self.nvalency
        else:
            return glob
        return gl




    # ----------------------------- PLOTTING -----------------------------


    def plot_mesh(self):
        setup_we_mpl()
        fig = plt.figure()
        ax  = fig.add_subplot(projection='3d')

        for el in self.elements:
            el.plot_element(ax)

        return fig, ax


    def plot_global_variable(self, var, fig=None, ax=None):
        setup_we_mpl()

        if ax==None:
            fig = plt.figure()
            ax  = fig.add_subplot(projection='3d')

        m = ax.scatter(self.gcoord_x, self.gcoord_y, self.gcoord_z, s=5, c=var)

        fig.colorbar(m)

        return ax