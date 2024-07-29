from abc import ABC
import numpy as np
from gll import gll
from hex import Hex

class Mesh():

    def __init__(self, nproc_xi, chunk=1, ngll=5):

        self.nproc_xi   = nproc_xi
        self.nproc_xip1 = nproc_xi + 1
        self.ngll       = ngll
        self.chunk      = chunk

        self.lats = None
        self.lons = None
        self.elements   = []

    def add_element(self,e):
        # Point to the mesh hex:
        e.hex = self.hex
        self.elements.append(e)

    def setup_gll(self):
        # get GLL points and weights:
        n = self.ngll
        self.gll, self.wgll = gll(n-1)

        # Total number of independent GLL points
        self.npoints_xi  =  self.nproc_xip1 + self.nproc_xi * (n - 2)
        self.npoints     = (self.npoints_xi**2) * n
        self.nelmts      = self.nproc_xi**2
        self.ibool       = np.zeros((n, n, n, self.nelmts), int)

        iglob = np.arange(self.npoints) + 1
        self.iglob = iglob.reshape((self.npoints_xi, self.npoints_xi, n), order='F')  # only 1 element deep

        i1 = self.iglob[:10, :7, 0]
        i2 = self.iglob[:, :, 1]
        i3 = self.iglob[:, :, 2]
        i4 = self.iglob[:, :, 3]
        i5 = self.iglob[:, :, 4]

        ielmt = 0
        for i in range(1, self.nproc_xip1):
            for j in range(self.nproc_xi):
                for izeta in range(n):
                    for ieta in range(n):
                        for ixi in range(n):

                            iind = i*(n-1) -  ixi
                            jind = j*(n-1) + ieta
                            kind = izeta

                            aa = self.iglob[iind,  jind,  kind]
                            self.ibool[ixi, ieta, izeta, ielmt] = aa
                ielmt += 1


        for elem in self.elements:
            elem._set_gll(self.gll, self.wgll)


    def create_hex(self):
        self.hex = Hex(self.ngll)


