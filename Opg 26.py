import numpy as np
import time as t
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
        p_vec[k] = np.trapz(p_func(x_int)*phi_k,x_int)  # integrer p(x)*phi_k med np.trapezoid
    return p_vec

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

# opg a
    # her skal der beregnes det samme som b men med kode-funktion fra opg 16 (viktor)

# opg b
def a_func(x):
    return x+1
def p_func(x):
    return 1

print(f"{'N':>6}  {'Sinusbasis (s)':>16}  {'Formfkt. (s)':>14}")
print("-" * 42)
for N in [10, 50, 100, 200]:
    nodes = [ a*(1/(N+1)) for a in range(N+2) ]
    times = []
    for _ in range(100):
        start = t.perf_counter() # start time counter

        K = assemble_K_general(nodes, a_func)  # saml stivhedsmatrixen
        p_vec = assemble_p_general(nodes, p_func)  # saml lastvektoren
        u = np.linalg.solve(K, p_vec)

        end = t.perf_counter()   # end time counter
        times.append(end-start)

    print(f"Mean time for N={N}: ", np.mean(times))
    print(f"Std for N={N}: ", np.std(times))

# plot FEM-løsningen sammen med den eksakte løsning?

# opg c
    # check plots

# opg d
    # ud fra c samt sammenlign a og b