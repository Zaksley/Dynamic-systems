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
def draw(x, y, t0, tf, modele):
    plt.plot(x, y, color='b')
    plt.xlabel("Temps")
    plt.ylabel("Population")
    plt.title("Dynamique des populations - Modèle " + modele)
    plt.show()

def malthus_simulation(y0, b, d, t0, tf):
    
    # Method
    meth = step_runge_kutta

    # Initialisation 
    epsilon = 10**-10
    
    def f(y, t):
        return (b-d)*y

    resolve = meth_epsilon(y0, t0, tf, epsilon, f, meth)
    draw(resolve[1], resolve[0], t0, tf, "Malthus")

    return resolve

def verhulst_simulation(y0, croissance, max, t0, tf):

    # Method
    meth = step_runge_kutta

    # Initialisation 
    epsilon = 10**-10
    
    def f(y, t):
        return y*croissance*(1 - y/max)

    resolve = meth_epsilon(y0, t0, tf, epsilon, f, meth)
    draw(resolve[1], resolve[0], t0, tf, "Verhulst")

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

def draw_2D(x, y_sheep, y_wolf):
    plt.plot(x, y_sheep, color='b')
    plt.plot(x, y_wolf, color='orange')
    plt.legend(["Proies", "Predateurs"])
    plt.xlabel("Temps")
    plt.ylabel("Population")
    plt.title("Dynamique des populations - Modèle Lotka-Volterra")
    plt.show() 

def draw_PP(y_sheep, y_wolf, y_sheep2, y_wolf2):
    plt.plot(y_sheep, y_wolf, color='b')
    plt.plot(y_sheep2, y_wolf2, color='r')
    plt.legend(["Modèle 1", "Modèle 2"])
    plt.xlabel("Population des proies")
    plt.ylabel("Population des prédateurs")
    plt.title("Variation du couple (proies, prédateurs)")
    plt.show() 

def draw_multiple_PP(res, size):
    for i in range(0, size-1, 2):
        if (len(res[i]) == len(res[i+1])):
            plt.plot(res[i], res[i+1])

    plt.xlabel("Population des proies")
    plt.ylabel("Population des prédateurs")
    plt.title("Variation du couple (proies, prédateurs)")
    plt.show() 


def lokta_simulation(y0, reprod_sheep, death_sheep, reprod_wolf, death_wolf, t0, tf):
      # Method
    meth = step_runge_kutta

    # Initialisation 
    epsilon = 10**-10
    
    def f(y, t):
        return np.array( [y[0] * (reprod_sheep - y[1]*death_sheep), y[1] * (reprod_wolf*y[0] - death_wolf) ]  )
 
    (resolve_y, resolve_x, norm) = meth_epsilon(y0, t0, tf, epsilon, f, meth)

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
                count += 1
        k += 1
    print(t2)
    print(t1)
    return t2- t1


def draw_multiple_simulation(y0, a, b, c, d, pas, i, t0, tf):
    
    colors_sheep = ["greenyellow", "lawngreen", "chartreuse", "springgreen", "lime" , "limegreen", "forestgreen", "green", "darkgreen"]
    colors_wolf = ['lightcoral', "coral", "salmon", "indianred", "orangered", "brown", "firebrick", "maroon", "red", "darkred"]

    res = []

    for j in range(i):
        y0 += pas
        (resolve_y, resolve_x) = lokta_simulation(y0, a, b, c, d, t0, tf)
        
        res.append([resolve_y[k][0] for k in range(len(resolve_y))])
        res.append([resolve_y[k][1] for k in range(len(resolve_y))])

        plt.plot(resolve_x,[resolve_y[k][0] for k in range(len(resolve_y))], color=colors_sheep[j])
        plt.plot(resolve_x, [resolve_y[k][1] for k in range(len(resolve_y))], color=colors_wolf[j])
    
    plt.legend(["Proies", "Predateurs"])
    plt.xlabel("Temps")
    plt.ylabel("Population")
    plt.title("Dynamique des populations - Modèle Lotka-Volterra")
    plt.show()

    draw_multiple_PP(res, i*2)


# Points singuliers
#
# (0,0) + (d/c, a/b)
#


if __name__ == "__main__":

    t0 = 0
    tf = 15

    # === Malthus 
    y0 = 2
    b = 2
    d = 1.3

    #malthus_simulation(y0, b, d, t0, tf)

    # === Verhulst
    croissance = b - d
    max = 2000

    #verhulst_simulation(y0, croissance, max, t0, tf)


    # === Interesting values (Lotka-Volterra) === 


    a = 2/3
    b = 4/3
    c = 1
    d = 1
    y0 = np.array([1.0, 1.0])
    t0 = 0
    tf = 45

            #Lotka simulation
    #(resolve_y, resolve_x) = lokta_simulation(y0, a, b, c, d, t0, tf)
    #draw_2D(resolve_x, [resolve_y[k][0] for k in range(len(resolve_y))], [resolve_y[k][1] for k in range(len(resolve_y))]); 

    a = 2/4
    b = 2/4
    c = 1/2
    d = 1/4
    y0 = np.array([2.0, 1.0])
    t0 = 0
    tf = 60

    #(resolve_y2, resolve_x2) = lokta_simulation(y0, a, b, c, d, t0, tf)

    #draw_PP([resolve_y[k][0] for k in range(len(resolve_y))], [resolve_y[k][1] for k in range(len(resolve_y))], 
    #       [resolve_y2[k][0] for k in range(len(resolve_y2))], [resolve_y2[k][1] for k in range(len(resolve_y2))])


        #Period
    #period_sheep = period([resolve_y[k][1] for k in range(len(resolve_y))], resolve_x)
    #print(period_sheep)

    # Points singuliers
    a = 2/4
    b = 2/4
    c = 1/2
    d = 1/4
    y0 = np.array([d/c, a/b])
    t0 = 0
    tf = 45

    pas = 0.05
    i = 8
    draw_multiple_simulation(y0, a, b, c, d, pas, i, t0, tf)


    
