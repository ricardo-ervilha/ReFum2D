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
x_vals = np.array([460, 942, 1992, 3962, 7826, 15644], dtype=float)
y_vals_mat = np.array([
    [ 0.0750389,
    0.027179,
    0.0131993,
    0.0053091,
    0.00408087,
    0.00219086
    ],
    [
    0.029763,
    0.011186,
    0.00677152,
    0.00221416,
    0.00370748,
    0.0017864
],
    [0.208114,
    0.111619,
    0.0661406,
    0.0383191,
    0.0412647,
    0.0223628
    ]
])

cmap = plt.get_cmap("tab10")
y_vals = np.sqrt(y_vals_mat[0]**2+ y_vals_mat[1]**2)
# Transformação: 1 / sqrt(x)
x_transformed = np.sqrt(4/x_vals)

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
plt.loglog(x_transformed, y_vals, marker='*')
plt.loglog(10**logx_fit, 10**logy_fit, label=label, linestyle='--', color = cmap(0.3))

plt.xlabel("$h$", fontsize=12)
plt.ylabel("Erro", fontsize=12)
plt.grid(True, which="both", linestyle="--", alpha=0.6)
plt.legend()
plt.savefig("./velocity_convergence.pdf", format='pdf', dpi=300)
plt.tight_layout()
plt.close()

cmap = plt.get_cmap("tab10")
y_vals = y_vals_mat[2]
# Transformação: 1 / sqrt(x)
x_transformed = np.sqrt(4/x_vals)

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
plt.loglog(x_transformed, y_vals, marker='*')
plt.loglog(10**logx_fit, 10**logy_fit, label=label, linestyle='--', color = 'orange')

plt.xlabel("$h$", fontsize=12)
plt.ylabel("Erro", fontsize=12)
plt.grid(True, which="both", linestyle="--", alpha=0.6)
plt.legend()
plt.savefig("./pressure_convergence.pdf", format='pdf', dpi=300)
plt.tight_layout()
plt.close()

