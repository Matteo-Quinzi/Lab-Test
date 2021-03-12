def Morse(x, D_e, a, x0):
    return D_e * (1. - exp(-a * (x - x0)))**2.