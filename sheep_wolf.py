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
def draw(array, t0, tf):
    n = (tf-t0)/(len(array))
    x = np.arange(t0, tf, n)

    plt.plot(x, array)
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
    draw(resolve, t0, tf)

    return resolve

print(malthus_simulation(300, 5, 4))

