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
        self.nelmts = 0

    def add_element(self,e):
        # Point to the mesh hex:
        e.hex = self.hex
        self.elements.append(e)
        self.nelmts += 1

    def setup_ibool(self):
        # get GLL points and weights:
        n = self.ngll
        self.gll, self.wgll = gll(n-1)
        self.range_ngll = range(self.ngll)

        # Total number of independent GLL points
        self.npoints_xi  =  self.nproc_xip1 + self.nproc_xi * (n - 2)
        self.npoints     = (self.npoints_xi**2) * n
        self.nelmts_hyp  = self.nproc_xi**2
        self.ibool       = np.zeros((n, n, n, self.nelmts_hyp), int)

        iglob = np.arange(self.npoints) + 1
        self.iglob = iglob.reshape((self.npoints_xi, self.npoints_xi, n), order='F')  # only 1 element deep

        ielmt = 0
        for i in range(self.nproc_xi):
            for j in range(self.nproc_xi):
                for izeta in range(n):
                    for ieta in range(n):
                        for ixi in range(n):

                            iind = i*(n-1) +  ixi
                            jind = j*(n-1) + ieta
                            kind = izeta

                            aa = self.iglob[iind,  jind,  kind]
                            self.ibool[ixi, ieta, izeta, ielmt] = aa
                ielmt += 1



    def create_hex(self):
        self.hex = Hex(self.ngll)


    def link_ibool_to_element(self):

        for iel in range(self.nelmts):
            self.elements[iel].ibool = self.ibool[:,:,:,iel]
