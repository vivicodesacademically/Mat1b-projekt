import numpy as np
"""
Immport funktioner herfra til benyttelse i forskellige opgaver.

En specifik funktion kan importes på følgende vis:
from FEM_functions import <function_name>

Alle funktioner kan importeres således:
from FEM_functions import *
"""


def assemble_K_sinus(N, a_func):
    """Stivhedsmatrix for sinusbasis med generel a(x)."""
    K = np.zeros((N, N))
    x = np.linspace(0, 1, 1000)


    for k in range(1, N + 1):
        dpsi_k = np.cos(k*np.pi*x)*k*np.pi  # den afledede af sin(k*pi*x)
        for l in range(1, N + 1):
            dpsi_l = np.cos(l*np.pi*x)*l*np.pi  # den afledede af sin(l*pi*x)
            try:
                K[k-1, l-1] = np.trapz(a_func(x)*dpsi_k*dpsi_l,x)  # integrer a(x)*dpsi_k*dpsi_l med np.trapezoid
            except:
                K[k-1, l-1] = np.trapezoid(a_func(x)*dpsi_k*dpsi_l,x)  # integrer a(x)*dpsi_k*dpsi_l med np.trapezoid
    return K

def assemble_p_sinus(N, p_func):
    """Lastvektor for sinusbasis med generel p(x)."""
    x = np.linspace(0, 1, 1000)
    p_vec = np.zeros(N)
    for k in range(1, N + 1):
        psi_k = np.sin(k*np.pi*x)  # beregn psi_k = sin(k*pi*x)
        try:
            p_vec[k-1] = np.trapz(p_func(x)*psi_k,x)  # integrer p(x)*psi_k med np.trapezoid
        except:
            p_vec[k-1] = np.trapezoid(p_func(x)*psi_k,x)  # integrer p(x)*psi_k med np.trapezoid
    return p_vec

### Opgave 25 ###
def assemble_p_general(nodes, p_func=None):
    """Saml lastvektoren p for vilkårligt net og kildefunktion p(x)."""
    N = len(nodes) - 2
    p_vec = np.zeros(N)
    if p_func is None:
        p_func = lambda x: 1.0
    for k in range(N):
        x_left = nodes[k]  # venstre endepunkt af phi_k's support
        x_right = nodes[k+2]  # højre endepunkt af phi_k's support
        x_int = np.linspace(x_left, x_right, 200)
        x_k = nodes[k+1]
        h_left = x_k-x_left  # bredde af venstre delelement
        h_right = x_right-x_k  # bredde af højre delelement
        phi_k = np.where(
            x_int <= x_k,
            (x_int - x_left) / h_left,
            (x_right - x_int) / h_right
        )
        try:
            p_vec[k] = np.trapz(p_func(x_int)*phi_k,x_int)  # integrer p(x)*phi_k med np.trapezoid
        except:
            p_vec[k] = np.trapezoid(p_func(x_int)*phi_k,x_int)  # integrer p(x)*phi_k med np.trapezoid
        
    return p_vec

### Opgave 24 ###
def assemble_K_general(nodes, a_func=None):
    """Saml K for vilkårligt net og funktion a(x)."""
    N = len(nodes) - 2
    K = np.zeros((N, N))
    if a_func is None:
        a_func = lambda x: 1.0
    for e in range(N+1):
        h_e = nodes[e+1] - nodes[e]  # elementbredde
        x_mid = 0.5*(nodes[e+1] + nodes[e])  # midtpunkt af elementet
        a_e = a_func(x_mid)  # a(x) evalueret i midtpunktet
        i = e - 1
        j = e
        # opdater K[i,i], K[j,j], K[i,j], K[j,i]
        # Husk at tjekke om i og j er indre knudepunkter
        k_local = (a_e / h_e) * np.array([[1, -1],
                                  [-1, 1]])
        global_nodes = [e, e+1]
        for a in range(2):
            for b in range(2):
                i_global = global_nodes[a]
                j_global = global_nodes[b]
                i_local = i_global - 1
                j_local = j_global - 1
                if 0 <= i_local < N and 0 <= j_local < N:
                    K[i_local, j_local] += k_local[a, b]
    return K

def FEM_y(x_plot, u):
    """Calculate the function values with formfunctions based on the u from the solved system."""
    N = len(u)
    h = 1.0 / (N + 1)  # skridtlængde
    phi = [0]*N 
    for k in range(1, N + 1):
        x_k = k*h  # knudepunktets position
        x_pk = (k+1)*h
        x_mk = (k-1)*h
        phi[k-1] = u[k-1]*np.piecewise(x_plot, 
            [(x_plot >= x_mk) & (x_plot < x_k), (x_plot < x_pk) & (x_plot >= x_k)], 
            [lambda x: ((x-x_mk)/(x_k-x_mk)), lambda x: ((x_pk-x)/(x_pk-x_k)), 0]) 
    return sum(phi)

def SinusBasis_y(x,u):
    """Calculate the function values with sinusbasis based on u from the solved system."""
    u_hat = np.zeros_like(x)

    for k in range(1, len(u) + 1):
        u_hat += u[k-1] * np.sin(k * np.pi * x)
    return u_hat