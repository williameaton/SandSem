from abc import ABC

class Mesh():

    def __init__(self):

        self.lats = None
        self.lons = None
        self.elements   = []

    def add_element(self,e):
        self.elements.append(e)

class Element():
    def __init__(self):
        self.cnodes_xyz = []
        self.cnodes     = []


    def add_corner_node(self,node):
        self.cnodes.append(node)





class Node():
    def __init__(self):
        self.global_id   = None
        self.element_id  = None
        self.alpha       = None
        self.beta        = None
        self.gamma       = None
        self.corner      = False