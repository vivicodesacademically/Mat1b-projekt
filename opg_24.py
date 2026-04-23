def assemble_K_general(nodes, a_func=None):
    """Saml K for vilkårligt net og funktion a(x)."""
    import numpy as np
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

N = 5
nodes = [ a*(1/(N+1)) for a in range(N+2) ]
print(assemble_K_general(nodes))