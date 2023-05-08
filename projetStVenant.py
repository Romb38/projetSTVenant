# -*- coding: utf-8 -*-

"""
Notes :
    U: [h,q]
    On considère que la hauteur d'eau et le débit sont constant aux bords
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sps
#import copy as cp

"""
Constantes informatiques

Le nombre de maille de notre projet correspond a LARGEUR * TPS_FIN
"""

LARGEUR = 100 #Valeur test
TPS_FINAL = 1000 #Valeur test
HAUTEUR_INITIALE = 2 #Valeur test


"""
Constantes physiques
"""
g = 9.81 #Constante de gravitée

"""
Matrice des hauteurs / débits
Les x_i sont les colones
Les t_n sont les lignes
"""

#mat_h = [[0 for i in range(LARGEUR)] for i in range (TPS_FINAL)]
mat_h = []
mat_q = []

"""
Définition de la grille d'étude
"""
#x_i = [ (i/LARGEUR) for i in range(0,NB_INTERVAL_POS)]
#t_n = [ (i/TPS_FINAL) for i in range(0,NB_INTERVAL_TPS)]

x_i = [ i for i in range(0,LARGEUR)]
t_n = [ (i/TPS_FINAL) for i in range(0,TPS_FINAL)]

delta_x = 1
delta_t = 0.1


"""
Condition CFL
"""
CFL = [0]

"""
Mise en place des condition initiales
"""

def loi_normale(mu,sgm,x):
    """
    Loi normale, pour la création de vagues
    """
    coeff1 = 1/(sgm*np.sqrt(2*np.pi))
    coeff2 = -1/2 * ((x-mu)/sgm)**2
    return coeff1 * np.exp(coeff2)


def h_initial(x):
    """
    Fonction de la hauteur au temps t=0
    
    Fonction test
    """
    return loi_normale(50,4,x) + HAUTEUR_INITIALE


def u_inital(x):
    """
    Fonction de la vitesse d'écoulement de l'eau au temps t=0
    
    Fonction test
    """
    return 0.5

def etat_initial_bords():
    """
    Mets la hauteurs au bords constante
    """
    for n in range(0,TPS_FINAL):
        mat_h[n][0] = HAUTEUR_INITIALE
        mat_h[n][LARGEUR-1] = HAUTEUR_INITIALE
    return
        
def ajout_ligne():
    """
    Ajoute une ligne a la matrice des hauteurs
    """
    mat_h.append([0 for i in range(LARGEUR)])
    mat_h[-1][0] = HAUTEUR_INITIALE
    mat_h[-1][LARGEUR-1] = HAUTEUR_INITIALE

    mat_q.append([1 for i in range(LARGEUR)])
    return

"""
Fonctions de calculs
"""

def F(U):
    """
    Calcul du flux
    Formule en dessous de la formule n°2 du papier de Maria
    """
    h = U[0]
    q = U[1]
    return [q, (q*q)/h + (g*h*h)/2]

def solveur_Rusanov(Ug,Ud):
    """
    Formule n°5 du papier de Maria
    """
    
    
    #Calcul du c selon la formule 6 du papier de Maria

    
    u_g = Ug[1]/Ug[0]
    h_g = Ug[0]
    
    u_d = Ud[1]/Ud[0]
    h_d = Ud[0]
    
    c_1 = abs(u_g) + np.sqrt(g*h_g)
    c_2 = abs(u_d) + np.sqrt(g*h_d)
    
    c = max(c_1,c_2)
    CFL[0] = c
    
    #Calcul du résultat voulu
    a_1 = (F(Ug)[0] + F(Ud)[0])/2 - c*(Ud[0]-Ug[0])/2
    a_2 = (F(Ug)[1] + F(Ud)[1])/2 - c*(Ud[1]-Ug[1])/2
    return [a_1,a_2]


def U_n(i,n):
    """
    Calcul de U a la position i au temps n+1
    Formule n°3 du papier de Maria
    """


    U_i1 = [mat_h[n][i-1],mat_q[n][i-1]]
    U_i2 = [mat_h[n][i+1],mat_q[n][i+1]]
    U_i = [mat_h[n][i],mat_q[n][i]]
    
    F_1 = solveur_Rusanov(U_i1,U_i)
    F_2 = solveur_Rusanov(U_i,U_i2)
    
    h_final = mat_h[n][i] - (delta_t/delta_x) * (F_2[0] - F_1[0])
    q_final = mat_q[n][i] - (delta_t/delta_x) * (F_2[1] - F_1[1])
    
    mat_h[n+1][i] = h_final
    mat_q[n+1][i] = q_final
    
    return
    
def main():
    
    ajout_ligne()

    #Mise en place des bords constants
    #etat_initial_bords()
    
    #Initialisation de la hauteur
    for i in range(1,LARGEUR-1):
        mat_h[0][i] = h_initial(i)
        mat_q[0][i] = h_initial(i) * u_inital(i)
    
    #Calcul des U_n pour toutes les positions/tous les temps
    t = 0
    n = 0
    while t < TPS_FINAL:
        t += CFL[0]
        ajout_ligne()
        for i in range(1,LARGEUR-1):
            U_n(i,n)
        n += 1
    
    #Affichage de l'état a la position finale
    
    for k in range(0,n,10):
        plt.plot(x_i,mat_h[k])
    plt.show()
    return 


def print_matrice(M):
    """Affiche une matrice correctement, FONCTION DE DEBUG"""
    for i in (M):
        print(M)


if __name__ == "__main__":
    main()