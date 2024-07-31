from abc import ABC
import numpy as np
from src.solver.timescheme import setup_timescheme

class Solver(ABC):
    def __init__(self, mesh, **kwargs):

        self.mesh = mesh

        # Default qwargs
        defaultKwargs = {'timescheme': None,
                         'print_kwargs': False}

        # Append to entered
        self.kwargs = {**defaultKwargs, **kwargs}

        if self.kwargs['print_kwargs']:
            self._print_kwargs()

        # Timescheme
        self.timescheme = self.kwargs['timescheme']
        if self.timescheme != None:
            self.ts = setup_timescheme(self, self.timescheme)



    def generate_variable_arrays(self, arrs, ndim=1):
        # Note: I believe in some cases the use of setattr is considered bad practice
        #       for security reasons. Here this seems like a low likelihood issue
        #       and it avoids the use of a dictionary, but maybe convert to a dictionary
        #       in the future
        if type(arrs) == str:
            arrs = [arrs]

        lena = len(arrs)

        if type(ndim)==int:
            # convert to a list:
            ndim = [ndim] * lena

        for a in range(lena):
            setattr(self, arrs[a], TransientVar(ndim=ndim[a], npoints=self.mesh.npoints))


    def zero_displacement(self,):
        # Zero all arrays:
        self.disp.prev *= 0
        self.disp.now  *= 0
        self.disp.next *= 0


    def generate_static_var(self, variable, ndim):
        # Note: I believe in some cases the use of setattr is considered bad practice
        #       for security reasons. Here this seems like a low likelihood issue
        #       and it avoids the use of a dictionary, but maybe convert to a dictionary
        #       in the future
        if ndim==1:
            vec = np.zeros((self.mesh.npoints))
        else:
            vec = np.zeros((ndim, self.mesh.npoints))
        setattr(self, variable, vec)


    def compute_spatial_scalar_gradient(self, sc):
        # Computes gradient of a global scalar field
        # sc is the scalar:
        # Local (element-wise) version of the scalar field
        #loc_sc = self.mesh.map_global_to_local(sc)

        m = self.mesh
        n = m.ngll

        # ldash ordering: col = gll nodes
        #                 row = order
        ldash = m.ldash

        # Local gradient
        grad = np.zeros((3, n, n, n, m.nelmts))

        for ielmt in range(m.nelmts):
            e    = m.elements[ielmt]
            ib   = e.ibool
            ijac = e.jacinv              # Jinv_{ij} is del Xi_i/del x_j
            for s in m.range_ngll:
                for t in m.range_ngll:
                    for n in m.range_ngll:
                        for idir in range(3):

                            for p in m.range_ngll:
                                grad[idir,s,t,n,ielmt] += sc[ib[p, t, n]] * ijac[0, idir, s, t, n] * ldash[p, s] + \
                                                          sc[ib[s, p, n]] * ijac[1, idir, s, t, n] * ldash[p, t] + \
                                                          sc[ib[s, t, p]] * ijac[2, idir, s, t, n] * ldash[p, n]

        grad_glob = np.zeros((3, m.npoints))
        for dim in range(3):
            grad_glob[dim,:] = self.mesh.map_local_to_global(local=grad[dim,:,:,:,:])

        return grad_glob



    def _print_kwargs(self):
        print('Kwargs: ', self.kwargs)




class TransientVar():
    def __init__(self, ndim, npoints):
        # In general how many different time points we need to store may vary
        # e.g. for a 5 point FD you would need to store the two previous timesteps
        # where as for an FD3 you would only need the previous 1
        self.prev = np.zeros((ndim, npoints))
        self.now  = np.zeros((ndim, npoints))
        self.next = np.zeros((ndim, npoints))