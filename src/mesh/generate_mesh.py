from src.mesh.mesh import Mesh

def generate_mesh(nchunks=1, **kwargs):
    defaultKwargs = {'ngll':      5,
                     'nproc_xi':  8,
                     'nelem_rad': 5,
                     'chunk':     1,
                     'rmin':    0.5,
                     'rmax':    1,
                     }

    kwargs = {**defaultKwargs, **kwargs}
    ngll = kwargs['ngll']



    if nchunks > 1:
        raise ValueError('Only implemented for one chunk currently')

    print(kwargs['nelem_rad'])

    m = Mesh(nproc_xi=kwargs['nproc_xi'], nelem_rad=kwargs['nelem_rad'], chunk=kwargs['chunk'], ngll=ngll)

    m.rmin = kwargs['rmin']
    m.rmax = kwargs['rmax']


    m.create_hex()
    m.setup_ibool()
    m.setup_elements()
    m.link_ibool_to_element()
    m.set_homogenous_property('cijkl', 3e10)
    m.set_homogenous_property('density', 2500)
    m.setup_integration()

    return m


