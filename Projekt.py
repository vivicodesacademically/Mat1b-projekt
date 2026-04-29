import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

import time

def assemble_K_sinus(N, a_func):
    """Stivhedsmatrix for sinusbasis med generel a(x)."""
    K = np.zeros((N, N))
    x = np.linspace(0, 1, 1000)


    for k in range(1, N + 1):
        dpsi_k = np.cos(k*np.pi*x)*k*np.pi  # den afledede af sin(k*pi*x)
        for l in range(1, N + 1):
            dpsi_l = np.cos(l*np.pi*x)*l*np.pi  # den afledede af sin(l*pi*x)
            K[k-1, l-1] = np.trapz(a_func(x)*dpsi_k*dpsi_l,x)  # integrer a(x)*dpsi_k*dpsi_l med np.trapezoid
    return K

def assemble_p_sinus(N, p_func):
    """Lastvektor for sinusbasis med generel p(x)."""
    x = np.linspace(0, 1, 1000)
    p_vec = np.zeros(N)
    for k in range(1, N + 1):
        psi_k = np.sin(k*np.pi*x)  # beregn psi_k = sin(k*pi*x)
        p_vec[k-1] = np.trapz(p_func(x)*psi_k,x)  # integrer p(x)*psi_k med np.trapezoid
    return p_vec

a_func = lambda x: np.where(x < 0.5, 2.0, 1.0)  # definer a(x) med np.where
p_func = lambda x: np.where(x < 3/5, 20, 0)  # definer p(x) med np.where

# (a) Saml K og p for N = 10
N = 10
K_sin = assemble_K_sinus(N, a_func)
p_sin = assemble_p_sinus(N, p_func)


# undersøg K, f.eks. med plt.spy eller print
print(K_sin)
print(p_sin)
# (b) Løs systemet og plot løsningen


# løs K*u = p og plot u_hat(x) = sum u_k*sin(k*pi*x)
u=np.linalg.solve(K_sin, p_sin)
print(np.vstack(u))

x_plot = np.linspace(0, 1, 200)
u_hat = np.zeros_like(x_plot)

for k in range(1, N + 1):
    u_hat += u[k-1] * np.sin(k * np.pi * x_plot)

print("u hat", u_hat)

plt.figure()
plt.plot(x_plot, u_hat, label=r' $\hat{u}(x)$')
# Plot den eksakte løsning u = 0.5*x*(1-x) hvis a=1, p=1
plt.xlabel('x')
plt.ylabel(r'$\hat{u}(x)$')
plt.legend()
plt.grid(True)
plt.title(f"Approksimative løsning, (N={N})")
plt.show()
"""u_hat = 0
for u_k in u:
    u_hat. lambda x: u[u_k]*np.sin(u_k*np.pi*x)
print(u_hat)"""