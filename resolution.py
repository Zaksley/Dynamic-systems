## Modules

import numpy as np
import matplotlib.pyplot as plt
from math import *

## Euler Method

# Computes step n+1 with Euler method
#
# @param y value at rank n
# @param t parameter at rank n
# @param f differential equation function
# @param h step size
# @return value at rank n+1
#
def step_euler(y, t, f, h):
    return np.array(y + h*f(y, t))


# Computes step n+1 with Middle point method
#
# @param y value at rank n
# @param t parameter at rank n
# @param f differential equation function
# @param h step size
# @return value at rank n+1
#
def step_middle_point(y, t, f, h):
    return np.array(y + h*f(y, t+(h/2)))

# Computes step n+1 with Heun method
#
# @param y value at rank n
# @param t parameter at rank n
# @param f differential equation function
# @param h step size
# @return value at rank n+1
#
def step_heun(y, t, f, h):
    p1 = f(y, t)
    y2 = y + h*p1
    p2 = f(y2, t+h)
    return np.array(y + h*(p1+p2)/2)

# Computes step n+1 with Heun method
#
# @param y value at rank n
# @param t parameter at rank n
# @param f differential equation function
# @param h step size
# @return value at rank n+1
#
def step_runge_kutta(y, t, f, h):
    p1 = f(y, t)
    y2 = y + h*p1/2
    p2 = f(y2, t+h/2)
    y3 = y + h*p2/2
    p3 = f(y3, t+h/2)
    y4 = y + h*p3
    p4 = f(y4, t+h)
    return np.array(y + h*(p1+2*p2+2*p3+p4)/6)

# Computes the value after N steps using a given method
#
# @param y0 value at initial parameter value
# @param t0 initial parameter value
# @param N number of steps
# @param h step size
# @param f differential equation function
# @param meth used method
# @return value at rank N
#
def meth_n_step(y0, t0, N, h, f, meth):
    y = np.array([0. for k in range(N)])
    y[0] = y0
    t = t0

    for i in range(N-1):
        y[i+1] = meth(y[i], t, f, h)
        t += h

    return y

# Computes the value after N steps using a given method
#
# @param y0 value at initial parameter value
# @param t0 initial parameter value
# @param tf parameter value at the end
# @param eps epsilon value
# @param f differential equation function
# @param meth used method
# @return approximation of y with error eps
#
def meth_epsilon(y0, t0, tf, eps, f ,meth):
    N = 10
    h = (tf - t0)/N
    i = 0

    yN = meth_n_step(y0, t0, N, h, f, meth) #Computed before the first iteration
    y2N = np.array([0])

    while(np.linalg.norm(yN-y2N[::2], inf) > eps and i <= 10000):
        yN = np.copy(y2N)
        y2N = meth_n_step(y0, t0, 2*N, h/2, f, meth)
        N *= 2
        h /= 2
        i += 1

    #print(N)
    return yN

## Tests

# Equation of dim 1 (solution : exp(-x)
# y(0) = 1
# y'(t) = -y(t)
def f_exp(y, t):
    return -y

# Equation of dim 1 (solution: exp(arctan(x))
# y(0) = 1
# y'(t) = y(t)/(1+t²)
def f_arcexp(y, t):
    return y/(1+t**2)


# Plots error between exp(-x) and its approximation step by step
# with only step_euler
def test__step_euler_exp():
    y = 1
    for i in range(100):
        x = i/100
        plt.scatter(x, abs(y-exp(-x)), c='b')
        y = step_euler(y, x, f_exp, 0.01)
    plt.yscale("log")
    plt.show()

# Equation of dim 2 (solution : [-sin(x), cos(x)]
# y0 = [1,0]
# y'(t) = [-y2(t), y1(t)] avec y(t) = [y1(t), y2(t)]
def f_sec(y, t):
    return np.array([-y[1], y[0]])

# Plots error between [cos(x), sin(x)] and its approximation step by step
# with only step_euler
def test__step_euler_dim2():
    y = np.array([1,0])
    for i in range(100):
        x = i/100
        #plt.scatter(y[0], y[1], c='r')
        #plt.scatter(cos(x), sin(x), c='b')
        plt.scatter(np.linalg.norm(abs(y - np.array([cos(x), sin(x)])), inf), x, c='b')
        y = step_euler(y, x, f_sec, 0.01)
    plt.yscale("log")
    plt.show()


# Plots error between exp(-x) and its apprxomiation with error 10e-4
def test__meth_epsilon_euler():
    euler = meth_epsilon(1, 0, 2, 10e-4, f_exp, step_euler)
    e = [exp(-i*(2/(len(euler)-1))) for i in range(len(euler))]
    x = [k*2/len(euler) for k in range(len(e))]

    plt.plot(x, [abs(euler[k] - e[k]) for k in range(len(e))])
    plt.yscale("log")
    plt.show()
    return


# Plots error between exp(arctan(x)) and its apprxomiation with error 10e-4
# using middle point method
def test__middle_point():
    euler = meth_epsilon(1, 0, 2, 10e-4, f_arcexp, step_middle_point)
    e = [exp(-i*(2/(len(euler)-1))) for i in range(len(euler))]
    x = [k*2/len(euler) for k in range(len(e))]

    plt.plot(x, [abs(euler[k] - e[k]) for k in range(len(e))])
    plt.yscale("log")
    plt.show()
    return

# Plots error between exp(arctan(x)) and its apprxomiation with error 10e-4
# using heun method
def test__heun():
    list = meth_epsilon(1, 0, 2, 10e-4, f_arcexp, step_euler)
    e = [exp(atan(i*(2/(len(list)-1)))) for i in range(len(list))]
    x = [k*2/len(list) for k in range(len(e))]

    plt.plot(x, [abs(list[k] - e[k]) for k in range(len(e))])
    plt.yscale("log")
    plt.show()
    return

# Plots error between exp(arctan(x)) and its apprxomiation with error 10e-4
# using runge_kutta method
def test__runge_kutta():
    list = meth_epsilon(1, 0, 2, 10e-4, f_arcexp, step_runge_kutta)
    e = [exp(atan(i*(2/(len(list)-1)))) for i in range(len(list))]
    x = [k*2/len(list) for k in range(len(e))]

    plt.plot(x, [abs(list[k] - e[k]) for k in range(len(e))])
    plt.yscale("log")
    plt.show()
    return


## Main

if __name__ == "__main__":
    test__runge_kutta()







