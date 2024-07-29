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
        f  = self.hex.faces
        c  = self.cnodes_xyz


        ax3D.scatter(c[:,0], c[:,1], c[:,2], s=0.5, alpha=0.4, c='k')

        for iface in range(6):
                conn = [x - 1 for x in f[iface].nodes]              # python indexing
                conn.append(conn[0])                                # joines up last edge
                ax3D.plot(c[conn,0], c[conn,1], c[conn,2],'k',      # plot
                          alpha=0.5, linewidth=0.3)


    def _set_gll(self, g, w):
        # g = array of gll
        # w = array of weights
        ngll = len(g)
        rngll = range(ngll)
        self.ibool_elem   = np.zeros((ngll,ngll,ngll))
        self.ibool_global = np.zeros((ngll,ngll,ngll))

