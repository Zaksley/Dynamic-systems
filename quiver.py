import matplotlib.pyplot as plt
import numpy as np
from math import *

def function_test(x, y):
    return np.array([  -2*x*(y[0]**2)  ])

def function_test_1D(x, y):
    return np.array([ y[0]/(1+x**2) ])

def draw_tangente_field_1D(F):
    X = np.arange(-2, 2, 0.1)
    Y = np.arange(-2, 2, 0.1)

    U = np.zeros( (len(X), len(Y)) )
    V = np.zeros( (len(X), len(Y)) )

    # Creation U & V matrices (Vectors)
    for i in range(len(X)):
        for j in range(len(Y)):
            x = X[i]
            y = Y[j]

            value = F(x, np.array([y]) )[0]

            # Vecteur (1, value)
            # Vecteur z = (1, value)
            # z = z/norme(z)
            # U = 1/norme(z)

            U[j][i] = (1 / sqrt(1 + value**2 ))
            V[j][i] = U[j][i] * value

    print(U)
    print(V)

    fig, ax = plt.subplots()
    q = ax.quiver(X, Y, U, V)
    #ax.quiverkey(q, X=0.3, Y=1.1, U=5,
    #        label='Quiver key, length = 10', labelpos='E')

    plt.show()


def function_test_2d(x, y): 
    return np.array([ -y, x ])


def draw_tangente_field_2D(F):
    X = np.arange(-2, 2, 0.1)
    Y = np.arange(-2, 2, 0.1)

    U = np.zeros( (len(X), len(Y)) )
    V = np.zeros( (len(X), len(Y)) )

    # Creation U & V matrices (Vectors)
    for i in range(len(X)):
        for j in range(len(Y)):
            x = X[i]
            y = Y[j]

            value = F(x, y)
            norme = sqrt(value[0]**2 + value[1]**2)

            U[j][i] = ( value[0] / norme)
            V[j][i] = ( value[1] / norme)

    print(U)
    print(V)

    fig, ax = plt.subplots()
    q = ax.quiver(X, Y, U, V)
    #ax.quiverkey(q, X=0.3, Y=1.1, U=5,
    #        label='Quiver key, length = 10', labelpos='E')

    plt.show()


if __name__ == "__main__":
    draw_tangente_field_1D(function_test_1D)
    draw_tangente_field_2D(function_test_2d)