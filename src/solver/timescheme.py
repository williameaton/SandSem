from abc import ABC

import numpy as np


def setup_timescheme(solver, ts):
    tsdict = {'CENTREDFD' : CentredFD(solver),
              'NEWMARK'   : Newmark(solver),
              'SYMPLECTIC': Symplectic(solver),
              }
    return tsdict[ts.upper()]


class Timescheme(ABC):
    def __init__(self, solver):
        self.solver = solver

def CentredFD(Timescheme):
    def __init__(self, solver):
        super().__init__()

def Newmark(Timescheme):
    def __init__(self, solver):
        super().__init__()

def Symplectic(Timescheme):
    def __init__(self, solver):
        super().__init__()




class TransientVar():
    def __init__(self, ndim, npoints):
        # In general how many different time points we need to store may vary
        # e.g. for a 5 point FD you would need to store the two previous timesteps
        # where as for an FD3 you would only need the previous 1
        self.prev = np.zeros((ndim, npoints))
        self.now  = np.zeros((ndim, npoints))
        self.next = np.zeros((ndim, npoints))