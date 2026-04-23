import numpy as np
import matplotlib.pyplot as plt
import time
from FEM_functions import *


a_func = lambda x: 1 + x
p_func = lambda x: 1 + 0*x


print(f"{'N':>6}  {'Sinusbasis (s)':>16}  {'Formfkt. (s)':>14}")
print("-" * 42)
for N in [10, 50, 100, 200, 5]:
    # Sinusbasis
    t0 = time.time()
    K_sin = assemble_K_sinus(N, a_func)   # saml med assemble_K_sinus
    p_sin = assemble_p_sinus(N, p_func)   # saml med assemble_p_sinus
    u_sin = np.linalg.solve(K_sin, p_sin)   # løs systemet
    t_sin = time.time() - t0

    # Formfunktioner på uniformt net
    nodes = [ a*(1/(N+1)) for a in range(N+2) ]   # uniformt net med N indre knudepunkter
    t0 = time.time()
    K_fem = assemble_K_general(nodes, a_func)   # saml med assemble_K_general
    p_fem = assemble_p_general(nodes, p_func)   # saml med assemble_p_general
    u_fem = np.linalg.solve(K_fem, p_fem)   # løs systemet
    t_fem = time.time() - t0

    print(f"{N:>6}  {t_sin:>16.2e}  {t_fem:>14.2e}")

# plot FEM-løsningen sammen med den eksakte løsning?
print(K_sin)
print(K_fem)
# opg c
    # check plots

x_plot = np.linspace(0, 1, 500)

plt.figure(figsize=(8, 4))
plt.xlabel('$x$')
plt.ylabel('$u(x)$')
plt.title(f'FEM-løsning plottet mod eksakte løsning')


plt.plot(x_plot, FEM_y(x_plot, u_fem),'r-', label=f'$\\phi_k, N={N}$')
plt.plot(x_plot, SinusBasis_y(x_plot, u_sin),'b--', label=f'$\\phi_k, N={N}$')

plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
# opg d
    # ud fra c samt sammenlign a og b