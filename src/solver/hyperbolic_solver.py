from src.solver.solver import Solver

class HyperbolicSolver(Solver):
    def __init__(self, mesh, **kwargs):
        super().__init__(mesh, **kwargs)

    def compute_forces(self):
        self.compute_Fvector()
        self.compute_KU()

    def compute_Fvector(self):
        pass

    def compute_KU(self):
        pass

    def plot_displacement(self):
        pass