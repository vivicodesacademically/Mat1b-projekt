import numpy as np
import matplotlib.pyplot as plt
from opg_24 import assemble_K_general
from opg_25 import assemble_p_general

#plot function
def FEM_y(x, u):
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

def a_func(x):
    if 0 <= x <= 1/2:
        return 2
    if 1/2 < x <= 1:
        return 1
def p_func(x):
    return 1

N=255
nodes = [ a*(1/(N+1)) for a in range(N+2) ]

x_plot=np.linspace(0,1,500)

#stykkevis sammesatte a
K = assemble_K_general(nodes, a_func)  # saml stivhedsmatrixen
p_vec = assemble_p_general(nodes, p_func)  # saml lastvektoren
u_a = np.linalg.solve(K, p_vec)

#a=1
K = assemble_K_general(nodes, lambda x:1)  # saml stivhedsmatrixen
p_vec = assemble_p_general(nodes, p_func)  # saml lastvektoren
u_a1 = np.linalg.solve(K, p_vec)

#a=2
K = assemble_K_general(nodes, lambda x:2)  # saml stivhedsmatrixen
p_vec = assemble_p_general(nodes, p_func)  # saml lastvektoren
u_a2 = np.linalg.solve(K, p_vec)

#plot af de tre løsninger
plt.figure(figsize=(8, 4))
plt.xlabel('$x$')
plt.ylabel('$u(x)$')
plt.title(f'FEM-løsning af stykkevis a(x) plottet mod hhv. a(x)=1 og a(x)=2')

#label=f'$a(x)=\begin{cases}2, & 0 \le x \le \frac{1}{2} \\ 1, & \frac{1}{2} < x \le 1\end{cases}, N={N}$' 
#label = f"a(x)=2 (0<=x<=0.5), 1 (0.5<x<=1), N={N}"
label = "a(x)=\n2,   0 ≤ x ≤ 1/2\n1,   1/2 < x ≤ 1"   
plt.plot(x_plot, FEM_y(x_plot, u_a),'r-', label=label)
plt.plot(x_plot,FEM_y(x_plot,u_a1),'b--', label=f'$a(x)=1, N={N}$')
plt.plot(x_plot,FEM_y(x_plot,u_a2),'g--', label=f'$a(x)=2, N={N}$')

plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()