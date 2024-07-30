from src.solver.timescheme import setup_timescheme, TransientVar

class HyperbolicSolver():
    def __init__(self, mesh, **kwargs):
        defaultKwargs = {'timescheme': None}
        kwargs = {**defaultKwargs, **kwargs}

        self.mesh = mesh
        self.timescheme = kwargs['timescheme']

        if self.timescheme!=None:
            self.ts = setup_timescheme(self, self.timescheme)


    def generate_variable_arrays(self, arrs):
        # Generates desired variable arrays:
        if 'displacement' in arrs:
            self.disp = TransientVar(ndim=1, npoints=self.mesh.npoints)


    def zero_displacement(self,):
        # Zero all arrays:
        self.disp.prev *= 0
        self.disp.now  *= 0
        self.disp.next *= 0

    def compute_forces(self):
        self.compute_Fvector()
        self.compute_KU()

    def compute_Fvector(self):
        pass

    def compute_KU(self):
        pass


    #def plot_displacement(self):