from src.solver.solver import Solver



class StaticSolver(Solver):

    def __init__(self, mesh, **kwargs):
        super().__init__(mesh, **kwargs)

        self._check_for_timescheme()


    def _check_for_timescheme(self):
        if self.timescheme != None:
            raise Warning(f"Timescheme '{self.timescheme}' was parsed but using static solver. Parsed timescheme is ignored.")
        self.timescheme = None