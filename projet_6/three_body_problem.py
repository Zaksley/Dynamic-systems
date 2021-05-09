## Modules

import numpy as np
import matplotlib.pyplot as plt
from math import *

## Tools

# Creates a gravitational force vector
#
# @param A First body position
# @param B Second body position
# @param mA mass of the first body
# @param mB mass of the second body
# @param G gravitationnal constant
# @return The force vector F{b->a}
#
#def gravitational_force(A, B, mA, mB, G):
#    d = ((A[0] - B[0])**2 + (A[1] - B[1])**2)**(3/2)
#    k = G*mA*mB
#    return np.array([k*(B[0]-A[0])/d, k*(B[1]-A[1])/d])

## Two-body problem

"""
Dynamic equation on body B :
a(t) = F{a->b}(t) / mB    = -G * x(t) * mA / (x(t)²+y(t)²)**(3/2)

We need 4 parameters :
mA, mB, x0A, x0B, v0B
(baucase v0A = 0; G = 1)

y(t) = [x(t), v(t)]
y(0) = [x0B, v0B]
y'(t) = [v(t), a(t)] = [x'(t), v'(t)]
"""

# Global variables:
global mA, mB, G
mA = 1
mB = 0.01
mC = 0.00000001
G = 1

# Differential equation
# y(t) = [x(t), v(t)]
# y(0) = [x0B, v0B]
# y'(t) = [v(t), a(t)] = [x'(t), v'(t)]
def two_body_func(y, t):
    return np.array([y[1], [-G*y[0][k]*mA / ((y[0][0]**2 + y[0][1]**2)**(3/2)) for k in range(2)]])


def two_body_problem(x0A, x0B, v0B, tf):
    y0 = np.array([x0B, v0B])

    (y, t, norm) = meth_epsilon(y0, 0, tf, 10e-4, two_body_func, step_runge_kutta)

    x = [y[k][0] for k in range(len(y))]
    posy = [x[k][1] for k in range(len(x))]
    posx = [x[k][0] for k in range(len(x))]

    plt.scatter(0,0)
    plt.plot(posx, posy)
    plt.axis('scaled')
    plt.show()

    return


## Three-body problem

# Differential equation
# y(t) = [xB(t), xC(t), vB(t), vC(t)]
# y(0) = [x0B, x0C, v0B, v0C]
# y'(t) = [vB(t), vC(t), aB(t), aC(t)]
def three_body_func(y, t):
    xb, yb, xc, yc = y[0][0], y[0][1], y[1][0], y[1][1]

    accxB = -G*xb*mA / ((xb**2 + yb**2)**(3/2))
    accxB += -G*(xb-xc)*mC / (((xc-xb)**2 + (yc-yb)**2)**(3/2))
    accyB = -G*yb*mA / ((xb**2 + yb**2)**(3/2))
    accyB += -G*(yb-yc)*mC / (((xc-xb)**2 + (yc-yb)**2)**(3/2))

    accxC = -G*xc*mA / ((xc**2 + yc**2)**(3/2))
    accxC += -G*(xc-xb)*mB / (((xc-xb)**2 + (yc-yb)**2)**(3/2))
    accyC = -G*yc*mA / ((xc**2 + yc**2)**(3/2))
    accyC += -G*(yc-yb)*mB / (((xc-xb)**2 + (yc-yb)**2)**(3/2))

    return np.array([y[2], y[3], [accxB, accyB], [accxC, accyC]])

def three_body_problem(x0A, x0B, x0C, v0B, v0C, tf):
    y0 = np.array([x0B, x0C, v0B, v0C])

    (y, t, norm) = meth_epsilon(y0, 0, tf, 10e-4, three_body_func, step_runge_kutta)

    xB = np.array([y[k][0] for k in range(len(y))])
    xC = np.array([y[k][1] for k in range(len(y))])
    posyB = [xB[k][1] for k in range(len(xB))]
    posxB = [xB[k][0] for k in range(len(xB))]
    posyC = [xC[k][1] for k in range(len(xC))]
    posxC = [xC[k][0] for k in range(len(xC))]

    plt.scatter(0,0)
    plt.plot(posxB, posyB)
    plt.plot(posxC, posyC)
    plt.axis('scaled')
    plt.show()

    return


## Main

if __name__ == "__main__":

    # Two-body problem (mA = 10, mB = 1)
    #two_body_problem(np.array([0,0]), np.array([2,0]), np.array([0,1]), 5) # Circular orbit
    #two_body_problem(np.array([0,0]), np.array([2,0]), np.array([0,4]), 5) # Too much speed
    #two_body_problem(np.array([0,0]), np.array([1,0]), np.array([0,0.5]), 5) # Not enough speed = too much acceleration
    #two_body_problem(np.array([0,0]), np.array([1,0]), np.array([0,0]), 5) # Straight fall

    #Restrained three-body problem (mA = 1, mB = 0.01, mC = 0.0001)
    three_body_problem(np.array([0,0]), np.array([1,0]), np.array([1.01,0]), np.array([0,1]), np.array([-0.01,0.02]), 10) # Double orbital
    # Too much speed to orbit around earth, but not enough to quit solar system :
    #three_body_problem(np.array([0,0]), np.array([1,0]), np.array([1.1,0]), np.array([0,1]), np.array([-0.1,0.2]), 10)
    # Too much speed to orbit around sun :
    #three_body_problem(np.array([0,0]), np.array([1,0]), np.array([1.1,0]), np.array([0,1]), np.array([-1,2]), 4)



