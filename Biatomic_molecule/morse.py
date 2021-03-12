from numpy import exp

def Morse(x):
    a = 1.
    x0 = 0.
    D_e = 1.
    return D_e * (1. - exp(-a * (x - x0)))**2.

print(Morse(1.))