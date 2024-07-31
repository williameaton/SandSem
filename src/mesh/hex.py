import numpy as np
from src.mesh.node import Node
from src.mesh.face import Face

class Hex():
    def __init__(self, ngll):

        self.cnodes = []

        # Create the corner nodes:
        for rz in range(2):
            rzngll = (ngll-1)*ngll*ngll*rz

            # NODE 1: ID 1  /  NODE 5: ID 1 + (NGLL-1)*NGLL*NGLL
            n = Node()
            n.ngll = ngll
            n.id   = 1 +  rzngll
            n.cn_id   = 1 + rz*4
            self.cnodes.append(n)

            # NODE 2: ID NGLL  /  NODE 6: ID NGLL + (NGLL-1)*NGLL*NGLL
            n = Node()
            n.ngll = ngll
            n.id   = ngll  + rzngll
            n.cn_id = 2 + rz*4
            self.cnodes.append(n)

            # NODE 3: ID NGLL x NGLL /  NODE 7: ID NGLL x NGLL + (NGLL-1)*NGLL*NGLL
            n = Node()
            n.ngll = ngll
            n.id   = ngll**2 + rzngll
            n.cn_id = 3 + rz*4
            self.cnodes.append(n)

            # NODE 4: NGLL(NGLL -1) +1  /  NODE 7: ID NGLL x (NGLL -1) + (NGLL-1)*NGLL*NGLL
            n = Node()
            n.ngll = ngll
            n.id   = ngll*(ngll-1)+1 + rzngll
            n.cn_id = 4 + rz*4
            self.cnodes.append(n)


        # Node ID ordering:
        self.node_order = np.arange(ngll**3).reshape((ngll,ngll,ngll), order='F')


        # Set up faces:
        self.faces = []
        self.faces.append(Face(1,  1, [1,2,6,5], self.node_order[-1, :, :]))   # Face 1: + YZ
        self.faces.append(Face(2,  1, [2,3,7,6], self.node_order[:, -1, :]))   # Face 2: + XZ
        self.faces.append(Face(3, -1, [4,3,7,8], self.node_order[0, :, :]))    # Face 3: - YZ
        self.faces.append(Face(4, -1, [1,4,8,5], self.node_order[:, 0, :]))    # Face 4: - XZ
        self.faces.append(Face(5, -1, [1,2,3,4], self.node_order[:, :, 0]))    # Face 5: - XY
        self.faces.append(Face(6,  1, [5,6,7,8], self.node_order[:, :, -1]))   # Face 6: + XY

