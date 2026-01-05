import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "palatino",
    "axes.labelsize": 14,
    "axes.titlesize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12
})

# Carrega o arquivo (valores separados por vírgula)
# Cada linha tem 3 floats
data = np.loadtxt("../outputs/flow_over_cylinder/Flow Over a cylinder_Re=25.txt", delimiter=",")

# data.shape -> (n_timesteps, 3)
t = np.arange(data.shape[0])  # eixo temporal

# Plot de cada coluna como uma série temporal
plt.figure()

plt.plot(t[::50], np.log10(data[::50, 0]), label=r"$u$")
plt.plot(t[::50], np.log10(data[::50, 1]), label=r"$v$")
plt.plot(t[::50], np.log10(data[::50, 2]), label=r"$p$")

plt.xlabel("Iterações")
plt.ylabel(r"$||\phi^{(n)} - \phi^{(n-1)}||_{\infty}$")
plt.legend()
plt.grid(True)

plt.savefig("iterationsFOC25.pdf", format='pdf', dpi=300)
