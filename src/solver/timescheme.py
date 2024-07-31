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


