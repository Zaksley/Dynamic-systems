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

def lokta_simulation(reprod_sheep, death_sheep, reprod_wolf, death_wolf):
    



if __name__ == "__main__":
    print(malthus_simulation(2, 2, 1)[0])
    print(verhulst_simulation(2, 2, 3)[0])