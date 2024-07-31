
class Boundary():
    def __init__(self, name):
        self.name     = name
        self.elements = []
        self.face     = None
        self.nfaces   = 0