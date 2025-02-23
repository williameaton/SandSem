# Some codes from Lion Krichner's Github:
import numpy as np

def gll(N):
    """
    Returns GLL (Gauss Lobato Legendre module with collocation points and
    weights)
    """
    # Initialization of integration weights and collocation points
    # [xi, weights] =  gll(N)
    # Values taken from Diploma Thesis Bernhard Schuberth
    if N == 2:
        xi = [-1.0, 0.0, 1.0]
        weights = [0.33333333, 1.33333333, 0.33333333]
    elif N == 3:
        xi = [-1.0, -0.447213595499957, 0.447213595499957, 1.0]
        weights = [0.1666666667, 0.833333333, 0.833333333, 0.1666666666]
    elif N == 4:
        xi = [-1.0, -0.6546536707079772, 0.0, 0.6546536707079772, 1.0]
        weights = [0.1, 0.544444444, 0.711111111, 0.544444444, 0.1]
    elif N == 5:
        xi = [-1.0, -0.7650553239294647, -0.285231516480645, 0.285231516480645,
              0.7650553239294647, 1.0]
        weights = [0.0666666666666667,  0.3784749562978470,
                   0.5548583770354862, 0.5548583770354862, 0.3784749562978470,
                   0.0666666666666667]
    elif N == 6:
        xi = [-1.0, -0.8302238962785670, -0.4688487934707142, 0.0,
              0.4688487934707142, 0.8302238962785670, 1.0]
        weights = [0.0476190476190476, 0.2768260473615659, 0.4317453812098627,
                   0.4876190476190476, 0.4317453812098627, 0.2768260473615659,
                   0.0476190476190476]
    elif N == 7:
        xi = [-1.0, -0.8717401485096066, -0.5917001814331423,
              -0.2092992179024789, 0.2092992179024789, 0.5917001814331423,
              0.8717401485096066, 1.0]
        weights = [0.0357142857142857, 0.2107042271435061, 0.3411226924835044,
                   0.4124587946587038, 0.4124587946587038, 0.3411226924835044,
                   0.2107042271435061, 0.0357142857142857]
    elif N == 8:
        xi = [-1.0, -0.8997579954114602, -0.6771862795107377,
              -0.3631174638261782, 0.0, 0.3631174638261782,
              0.6771862795107377, 0.8997579954114602, 1.0]
        weights = [0.0277777777777778, 0.1654953615608055, 0.2745387125001617,
                   0.3464285109730463, 0.3715192743764172, 0.3464285109730463,
                   0.2745387125001617, 0.1654953615608055, 0.0277777777777778]
    elif N == 9:
        xi = [-1.0, -0.9195339081664589, -0.7387738651055050,
              -0.4779249498104445, -0.1652789576663870, 0.1652789576663870,
              0.4779249498104445, 0.7387738651055050, 0.9195339081664589, 1.0]
        weights = [0.0222222222222222, 0.1333059908510701, 0.2248893420631264,
                   0.2920426836796838, 0.3275397611838976, 0.3275397611838976,
                   0.2920426836796838, 0.2248893420631264, 0.1333059908510701,
                   0.0222222222222222]
    elif N == 10:
        xi = [-1.0, -0.9340014304080592, -0.7844834736631444,
              -0.5652353269962050, -0.2957581355869394, 0.0,
              0.2957581355869394, 0.5652353269962050, 0.7844834736631444,
              0.9340014304080592, 1.0]
        weights = [0.0181818181818182, 0.1096122732669949, 0.1871698817803052,
                   0.2480481042640284, 0.2868791247790080, 0.3002175954556907,
                   0.2868791247790080, 0.2480481042640284, 0.1871698817803052,
                   0.1096122732669949, 0.0181818181818182]
    elif N == 11:
        xi = [-1.0, -0.9448992722228822, -0.8192793216440067,
              -0.6328761530318606, -0.3995309409653489, -0.1365529328549276,
              0.1365529328549276, 0.3995309409653489, 0.6328761530318606,
              0.8192793216440067, 0.9448992722228822, 1.0]
        weights = [0.0151515151515152, 0.0916845174131962, 0.1579747055643701,
                   0.2125084177610211, 0.2512756031992013, 0.2714052409106962,
                   0.2714052409106962, 0.2512756031992013, 0.2125084177610211,
                   0.1579747055643701, 0.0916845174131962, 0.0151515151515152]
    elif N == 12:
        xi = [-1.0, -0.9533098466421639, -0.8463475646518723,
              -0.6861884690817575, -0.4829098210913362, -0.2492869301062400,
              0.0, 0.2492869301062400, 0.4829098210913362,
              0.6861884690817575, 0.8463475646518723, 0.9533098466421639,
              1.0]
        weights = [0.0128205128205128, 0.0778016867468189, 0.1349819266896083,
                   0.1836468652035501, 0.2207677935661101, 0.2440157903066763,
                   0.2519308493334467, 0.2440157903066763, 0.2207677935661101,
                   0.1836468652035501, 0.1349819266896083, 0.0778016867468189,
                   0.0128205128205128]
    else:
        raise NotImplementedError

    return np.array(xi), np.array(weights)


def lagrange_alt(N, i, x, xi):
    """
    Function to calculate  Lagrange polynomial for order N
    and polynomial i [0, N] at location x at given collacation points xi
    (not necessarily the GLL-points)
    """
    fac = 1
    for j in range(0, N+1):
        if j != i:
            fac = fac * ((x - xi[j]) / (xi[i] - xi[j]))
    return fac


def lagrange(N, i, x):
    """
    Function to calculate  Lagrange polynomial for order N and polynomial
    i[0, N] at location x.
    """

    [xi, weights] = gll(N)
    fac = 1
    for j in range(-1, N):
        if j != i:
            fac = fac * ((x - xi[j + 1]) / (xi[i + 1] - xi[j + 1]))
    return fac




def lagrange1st(N):
    """
    # Calculation of 1st derivatives of Lagrange polynomials
    # at GLL collocation points
    # out = legendre1st(N)
    # out is a matrix with columns -> GLL nodes
    #                        rows -> order
    """
    out = np.zeros([N+1, N+1])

    [xi, w] = gll(N)

    # initialize dij matrix (see Funaro 1993 or Diploma thesis Bernhard
    # Schuberth)
    d = np.zeros([N + 1, N + 1])

    for i in range(-1, N):
        for j in range(-1, N):
            if i != j:
                d[i + 1, j + 1] = legendre(N, xi[i + 1]) / \
                    legendre(N, xi[j + 1]) * 1.0 / (xi[i + 1] - xi[j + 1])
            if i == -1:
                if j == -1:
                    d[i + 1, j + 1] = -1.0 / 4.0 * N * (N + 1)
            if i == N-1:
                if j == N-1:
                    d[i + 1, j + 1] = 1.0 / 4.0 * N * (N + 1)

    # Calculate matrix with 1st derivatives of Lagrange polynomials
    for n in range(-1, N):
        for i in range(-1, N):
            sum = 0
            for j in range(-1, N):
                sum = sum + d[i + 1, j + 1] * lagrange(N, n, xi[j + 1])

            out[n + 1, i + 1] = sum
    return(out)



def legendre(N, x):
    """
    Returns the value of Legendre Polynomial P_N(x) at position x[-1, 1].
    """
    P = np.zeros(2 * N)

    if N == 0:
        P[0] = 1
    elif N == 1:
        P[1] = x
    else:
        P[0] = 1
        P[1] = x
    for i in range(2, N + 1):
        P[i] = (1.0 / float(i)) * ((2 * i - 1) * x * P[i - 1] - (i - 1) *
                                   P[i - 2])

    return(P[N])



def lagrange2nd(hdd, xi, NGL):
    # Credit L Sawade
    hdd[:,:] = 0.0
    # Loop over each part of matrix:
    for j in range(NGL+1):
        for r in range(NGL+1):

            # begin calc matrix for r, j
            for i in range(NGL+1):


                s1 = 0.0
                if i!=j:
                    for m in range(NGL+1):
                        p1 = 1.0
                        if m!=i and m!=j:
                            for n in range(NGL+1):
                                if n!=i and n!=j and n!=m:
                                    p1 = p1 * (xi[r] - xi[n])/(xi[j] - xi[n])
                            s1 += (1.0/(xi[j]-xi[m])) * p1

                    hdd[r, j] += (1.0/ (xi[j] - xi[i]))*s1
    return hdd


def mapping_deriv_3GLL_analytical(x,y,z):
    # x, y, z should be 3 x 3 x 3 matrices with values of X,Y,Z at GLL points
    # Jij = del x_i/del Xi_j
    ngll = 3
    xi, w  = gll(ngll-1)
    J = np.zeros((3,3,ngll,ngll,ngll))

    # columns -> GLL nodes
    # rows    -> order
    ldash = lagrange1st(ngll-1)


    c = [x,y,z]
    # Loop over GLL points to compute derivatives at
    for s in range(ngll):
        for t in range(ngll):
            for n in range(ngll):
                # Loop over x,y,z
                for i in range(3):
                    cc = c[i]
                    # Loop over xi, eta, zeta
                    for j in range(3):
                        val = 0

                        for p in range(ngll):
                            if  j == 0:
                                val += cc[p,t,n]   * ldash[p,s]
                            elif j == 1:
                                val += cc[s, p, n] * ldash[p,t]
                            elif j == 2:
                                val += cc[s, t, p] *  ldash[p,n]
                            else:
                                raise ValueError()


                        J[i,j,s,t,n] = val


