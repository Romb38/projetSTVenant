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

mat_h = [[0 for i in range(LARGEUR)] for i in range (TPS_FINAL)]
mat_q = [[1 for i in range(LARGEUR)] for i in range (TPS_FINAL)]

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
Mise en place des condition initiales
"""

def loi_normale(mu,sgm,x):
    """
    Loi normale, pour la création de vagues
    """
    coeff1 = 1/(sgm*np.sqrt(2*np.pi))
    coeff2 = -1/2 * ((x-mu)/sgm)**2
    return coeff1 * np.exp(coeff2)

def fond_constant(x):
    """
    Fonction de topographie CONSTANT AU COURS DU TEMPS


    Fonction test
    """
    return np.cos(1/20*x)/14 + 1.8

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

def solveur_Rusanov(Ug,Ud,zG,zD,side = "D"):
    """
    Formule n°5 du papier de Maria
    """


    z_star = max(zG,zD) #Formule 24 + du papier de Maria

    
    #Calcul du c selon la formule 6 du papier de Maria

    
    u_g = Ug[1]/Ug[0]
    h_g = Ug[0]
    
    u_d = Ud[1]/Ud[0]
    h_d = Ud[0]


    #Formule 24 du papier de Maria
    hg_star = np.abs(h_g + zG - z_star)
    hd_star = np.abs(h_d + zD - z_star)
    Ud_star = [hd_star, hd_star*u_d]
    Ug_star = [hg_star, hg_star*u_g]

    c_1 = abs(u_g) + np.sqrt(g*hg_star)
    c_2 = abs(u_d) + np.sqrt(g*hd_star)
    
    c = max(c_1,c_2)
    
    #Calcul du résultat voulu
    a_1 = (F(Ug_star)[0] + F(Ud_star)[0])/2 - c*(Ud_star[0]-Ug_star[0])/2
    a_2 = (F(Ug_star)[1] + F(Ud_star)[1])/2 - c*(Ud_star[1]-Ug_star[1])/2

    #Formule 23 du papier de Maria
    if side == "D":
        a_2 += g/(2*(h_d)**2 - (hd_star)**2)
    else :
        a_2 += g/(2*(h_g)**2 - (hg_star)**2)

    return [a_1,a_2]


def U_n(i,n):
    """
    Calcul de U a la position i au temps n+1
    Formule n°3 du papier de Maria
    """

    zG = fond_constant(i)
    zD = fond_constant(i+1)

    U_i1 = [mat_h[n][i-1],mat_q[n][i-1]]
    U_i2 = [mat_h[n][i+1],mat_q[n][i+1]]
    U_i = [mat_h[n][i],mat_q[n][i]]
    
    F_1 = solveur_Rusanov(U_i1,U_i,zG,zD,"D")
    F_2 = solveur_Rusanov(U_i,U_i2,zG,zD,"G")
    
    h_final = mat_h[n][i] - (delta_t/delta_x) * (F_2[0] - F_1[0])
    q_final = mat_q[n][i] - (delta_t/delta_x) * (F_2[1] - F_1[1])
    
    mat_h[n+1][i] = h_final
    mat_q[n+1][i] = q_final
    
    return
    
def main():
    
    #Mise en place des bords constants
    etat_initial_bords()
    
    #Initialisation de la hauteur
    for i in range(1,LARGEUR-1):
        mat_h[0][i] = h_initial(i)
        mat_q[0][i] = h_initial(i) * u_inital(i)
    
    #Calcul des U_n pour toutes les positions/tous les temps
    for n in range (0,TPS_FINAL-1):
        for i in range(1,LARGEUR-1):
            U_n(i,n)
    
    #Affichage de l'état a la position finale
    
    for n in range(0,TPS_FINAL,10):
        plt.plot(x_i,mat_h[n])

    fond = [fond_constant(i) for i in range (0,LARGEUR)]
    plt.plot(x_i,fond)
    plt.show()
    return 


def print_matrice(M):
    """Affiche une matrice correctement, FONCTION DE DEBUG"""
    for i in (M):
        print(M)


if __name__ == "__main__":
    main()