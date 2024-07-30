import numpy as np

class Element():
    def __init__(self):
        self.cnodes_xyz = []
        self.prop = {}

    def add_corner_node(self,node):
        self.cnodes.append(node)


    def cnodes_xyz_to_np(self):
        # Converts a number of cnodes_xyz (hopefully 8) to a np array)
        cn = np.zeros((8,3))
        for i in range(8):
            for j in range(3):
                cn[i,j] = self.cnodes_xyz[i][j]
        self.cnodes_xyz = cn

    def plot_element(self, ax3D):
        cn  = self.cnodes_xyz
        cx = self.coord_x
        cy = self.coord_y
        cz = self.coord_z

        # Plot all nodes
        #ax3D.scatter(cx.flatten(),
        #             cy.flatten(),
        #r             cz.flatten(), s=2, alpha=0.4, c='k')

        # Plot corner nodes
        #ax3D.scatter(cn[:,0], cn[:,1], cn[:,2], s=2, alpha=0.4, c='r')
        a = 0.2
        # Plot edges:
        for yy in [0,-1]:
            for xx in [0, -1]:
                ax3D.plot(cx[xx, yy, :],
                          cy[xx, yy, :],
                          cz[xx, yy, :], 'k', alpha=a)

        for zz in [0,-1]:
            for xx in [0, -1]:
                ax3D.plot(cx[xx, :, zz],
                          cy[xx, :, zz],
                          cz[xx, :, zz], 'k', alpha=a)

        for zz in [0, -1]:
            for yy in [0, -1]:
                ax3D.plot(cx[:, yy, zz],
                          cy[:, yy, zz],
                          cz[:, yy, zz], 'k', alpha=a)


        # Plot nodes and or IBOOL ids of all GLL
        for k in [-1]: #self.mesh.range_ngll:
            for j in self.mesh.range_ngll:
                for i in self.mesh.range_ngll:
                    """ax3D.scatter(cx[i, j, k],
                                 cy[i, j, k],
                                 cz[i, j, k], alpha=0.5, s=5, c='b')"""
                    #ax3D.text(cx[i, j, k], cy[i, j, k], cz[i, j, k], s=str(self.ibool[i,j,k]), fontsize=5)

    def _setup_gll(self):
        print()



    def compute_jacobian(self):
        self.jac = np.zeros((3,3, self.ngll, self.ngll, self.ngll))
        self.jac3D = np.zeros((self.ngll, self.ngll, self.ngll))

        # columns -> GLL nodes
        # rows    -> order
        c = [self.coord_x, self.coord_y, self.coord_z]

        # Loop over GLL points to compute derivatives at
        for s in range(self.ngll):
             for t in range(self.ngll):
                for n in range(self.ngll):
                    # Jij = del x_i/del Xi_j
                    J = np.zeros((3, 3))

                    # Loop over x,y,z
                    for i in range(3):
                        cc = c[i]
                        # Loop over xi, eta, zeta
                        for j in range(3):
                            val = 0

                            for p in range(self.ngll):
                                if j == 0:
                                    val += cc[p, t, n] * self.mesh.ldash[p, s]
                                elif j == 1:
                                    val += cc[s, p, n] * self.mesh.ldash[p, t]
                                elif j == 2:
                                    val += cc[s, t, p] * self.mesh.ldash[p, n]
                                else:
                                    raise ValueError()
                            J[i,j] += val

                    # Store the jacobian
                    self.jac[:,:,s,t,n]   = J
                    self.jac3D[s,t,n] = np.linalg.det(J)
