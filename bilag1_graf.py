import matplotlib.pyplot as plt
import numpy as np

def exact_u(x):
    return x*(1-x)/2

x = np.linspace(0,1,100)
plt.figure(figsize=(8, 4))
plt.xlabel("x")
plt.ylabel("$u(x)$")
plt.plot(x, exact_u(x), 'r-', label="Titel")
plt.scatter(0.5, 0.125, color=(0,1,1), alpha=1.0, zorder=5,label="max = (0.5, 0.125)")
plt.annotate("max = (0.5, 0.125)",
             (0.5, 0.125),
             xytext=(0, -20),
             textcoords='offset points',
             color="cyan",
             ha='center')
plt.show()