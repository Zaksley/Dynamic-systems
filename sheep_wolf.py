import matplotlib.pyplot as plt
import numpy as np
from resolution import *
from quiver import *

# === #
#
#   Question (1)
#           Modele Mathus : Birthds - Deaths : Croissance uniquement lié à ces paramètres
#           Modèle Verhust : 1 - Pop_actuelle/Max pop : Facteur limitant (ressources, etc..)
#               - Si pop faible => Monte forte
#               - Si pop forte => Monte peu
#
# === #

# 2) 
def draw(x, y, t0, tf):
    plt.plot(x, y)
    plt.show()

def malthus_simulation(y0, b, d):
    
    # Method
    meth = step_runge_kutta

    # Initialisation 
    t0 = 0
    tf = 3
    epsilon = 10**-10
    
    def f(y, t):
        return (b-d)*y

    resolve = meth_epsilon(y0, t0, tf, epsilon, f, meth)
    draw(resolve[1], resolve[0], t0, tf)

    return resolve

def verhulst_simulation(y0, croissance, max):

    # Method
    meth = step_runge_kutta

    # Initialisation 
    t0 = 0
    tf = 3
    epsilon = 10**-10
    
    def f(y, t):
        return croissance*(1 - y/max)

    resolve = meth_epsilon(y0, t0, tf, epsilon, f, meth)
    draw(resolve[1], resolve[0], t0, tf)

    return resolve




# === #
#
#   Question (3)
#           Evolution sheep = Population sheep * (a = Coeff reproduction - b = bouffe_sheep * Population wolf)
#           Modèle très similaire à Malthus
#
#           Evolution wolf = Population wolf * (c = reproduction en fonction des proies rencontrés * Population sheep - death wolf)
#               
#           a = reproduction sheep (constant)
#           b = mortalité en fonction du nombre de wolfs
#           c = reproduction en fonction du nombre de sheeps
#           d = mortalité wolf (constant)
# === #

def draw_2D(x, y_sheep, y_wolf, t0, tf):
    plt.plot(x, y_sheep, 'b')
    plt.plot(x, y_wolf, 'r')
    plt.show() 


def lokta_simulation(y0, reprod_sheep, death_sheep, reprod_wolf, death_wolf):
      # Method
    meth = step_runge_kutta

    # Initialisation 
    t0 = 0
    tf = 40
    epsilon = 10**-10
    
    def f(y, t):
        return np.array( [y[0] * (reprod_sheep - y[1]*death_sheep), y[1] * (reprod_wolf*y[0] - death_wolf) ]  )
 
    (resolve_y, resolve_x, norm) = meth_epsilon(y0, t0, tf, epsilon, f, meth)
    print(resolve_y)
    draw_2D(resolve_x, [resolve_y[k][0] for k in range(len(resolve_y))], [resolve_y[k][1] for k in range(len(resolve_y))], t0, tf); 

    return (resolve_y, resolve_x)

def period(y,t):
    count=0
    k=1
    t1=0
    t2=0
    while (count < 2) and (k<len(y)-1):
        if y[k-1]<=y[k] and y[k]>=y[k+1]:
            if count==0:
                t1=t[k]
                count += 1
            else:
                t2=t[k]
                count =+ 1
        k += 1
    print(t2)
    print(t1)
    return t2- t1


# Solutions
#
# 
#   
#
#
#

# Points singuliers
#
# (0,0) + (d/c, a/b)
#


if __name__ == "__main__":
    #print(malthus_simulation(2, 2, 1)[0])
    #print(verhulst_simulation(2, 2, 3)[0])

    a = 2/3
    b = 4/3
    c = 1
    d = 1
    y0 = np.array([1, 1])

    (resolve_y, resolve_x ) = lokta_simulation(y0, a, b, c, d)
    period_sheep = period([resolve_y[k][0] for k in range(len(resolve_y))], resolve_x)
    print(period_sheep)