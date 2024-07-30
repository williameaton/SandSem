import numpy as np


class Element():
    def __init__(self):
        self.cnodes_xyz = []


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
        ax3D.scatter(cn[:,0], cn[:,1], cn[:,2], s=2, alpha=0.4, c='r')

        # Plot edges:
        for yy in [0,-1]:
            for xx in [0, -1]:
                ax3D.plot(cx[xx, yy, :],
                          cy[xx, yy, :],
                          cz[xx, yy, :], 'k')

        for zz in [0,-1]:
            for xx in [0, -1]:
                ax3D.plot(cx[xx, :, zz],
                          cy[xx, :, zz],
                          cz[xx, :, zz], 'k')

        for zz in [0, -1]:
            for yy in [0, -1]:
                ax3D.plot(cx[:, yy, zz],
                          cy[:, yy, zz],
                          cz[:, yy, zz], 'k')


        # Plot nodes and or IBOOL ids of all GLL
        """for k in range(1): #self.mesh.range_ngll:
            for j in self.mesh.range_ngll:
                for i in self.mesh.range_ngll:
                    ax3D.scatter(cx[i, j, k],
                                 cy[i, j, k],
                                 cz[i, j, k], alpha=0.5, s=5, c='b')"""
                    #ax3D.text(cx[i, j, k], cy[i, j, k], cz[i, j, k], s=str(self.ibool[i,j,k]), fontsize=5)

    def _setup_gll(self):

        # Create the GLL coordinates
        print()


