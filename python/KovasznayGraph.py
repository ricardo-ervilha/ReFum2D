import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.path import Path

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino"],
    "text.latex.preamble": r"\usepackage{mathpazo}"
})

import matplotlib.pyplot as plt
import numpy as np

# Dados fornecidos
x_vals = np.array([460, 942, 1992, 3962, 7826], dtype=float)
y_vals_mat = np.array([
    [ 0.0372087,
    0.0135201,
    0.00653967,
    0.00261017,
    0.0021225],
    [
   0.0155147,
    0.00589536,
    0.00356925,
    0.00114755,
    0.00195394
],
    [0.110161,
    0.0598823,
    0.0354144,
    0.0209899,
    0.0224716]
])

cmap = plt.get_cmap("tab10")
y_vals = np.sqrt(y_vals_mat[0]**2 + y_vals_mat[1]**2)
# Transformação: 1 / sqrt(x)
x_transformed = 1 / np.sqrt(x_vals)

# Regressão linear em escala log-log
logx = np.log10(x_transformed)
logy = np.log10(y_vals)

# Ajuste linear
slope, intercept = np.polyfit(logx, logy, 1)

#   print("Coeficiente angular =", slope)

# Reta ajustada no espaço log-log
logx_fit = np.linspace(logx.min(), logx.max(), 200)
logy_fit = slope * logx_fit + intercept

order = fr"h^{{{slope:.2f}}}"
label = fr"$\vec{{u}}$, $\mathcal{{O}}({order})$"

# Gráfico
plt.figure(figsize=(8,6))
plt.loglog(x_transformed, y_vals, '*', markersize=8,  markeredgewidth=1.5)
plt.loglog(10**logx_fit, 10**logy_fit, label=label, color = cmap(0.3))

plt.xlabel("$h$", fontsize=12)
plt.ylabel("Erro", fontsize=12)
plt.grid(True, which="both", linestyle="--", alpha=0.6)
plt.legend()
plt.show()
plt.tight_layout()
plt.close()

cmap = plt.get_cmap("tab10")
y_vals = y_vals_mat[2]
# Transformação: 1 / sqrt(x)
x_transformed = 1 / np.sqrt(x_vals)

# Regressão linear em escala log-log
logx = np.log10(x_transformed)
logy = np.log10(y_vals)

# Ajuste linear
slope, intercept = np.polyfit(logx, logy, 1)

#   print("Coeficiente angular =", slope)

# Reta ajustada no espaço log-log
logx_fit = np.linspace(logx.min(), logx.max(), 200)
logy_fit = slope * logx_fit + intercept

order = fr"h^{{{slope:.2f}}}"
label = fr"$p$, $\mathcal{{O}}({order})$"

# Gráfico
plt.figure(figsize=(8,6))
plt.loglog(x_transformed, y_vals, '*', markersize=8,  markeredgewidth=1.5)
plt.loglog(10**logx_fit, 10**logy_fit, label=label, color ='orange')

plt.xlabel("$h$", fontsize=12)
plt.ylabel("Erro", fontsize=12)
plt.grid(True, which="both", linestyle="--", alpha=0.6)
plt.legend()
plt.show()
plt.tight_layout()
plt.close()

