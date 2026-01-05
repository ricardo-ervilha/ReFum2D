import meshio
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

# ======================================================
# 1. Ler o VTK
# ======================================================
mesh = meshio.read("../outputs/backward_facing_step/Backward Facing Step_Re=100.vtk")

points = mesh.points[:, :2]
cells = mesh.cells_dict["triangle"]

velocity_c = mesh.cell_data["velocity"][0]
u_c = velocity_c[:, 0]

# ======================================================
# 2. Converter cell data -> pseudo-point data
#    (necessário para interpolar nas arestas)
# ======================================================
n_points = points.shape[0]
u_p = np.zeros(n_points)
count = np.zeros(n_points)

for cid, cell in enumerate(cells):
    for nid in cell:
        u_p[nid] += u_c[cid]
        count[nid] += 1

u_p /= np.maximum(count, 1)

# ======================================================
# 3. Interseção das arestas com x = 1.08
# ======================================================
x_cut = 3.53
tol = 1e-12

x = np.array([
    3.841463414634146,
    3.841463414634146,
    4.390243902439023,
    4.390243902439023,
    5.121951219512194,
    5.853658536585366,
    7.134146341463413,
    8.414634146341461,
    10.243902439024389,
    12.073170731707314,
    13.536585365853655,
    14.634146341463412,
    15.365853658536581,
    15.914634146341461,
    16.097560975609756,
    14.634146341463412,
    13.902439024390242,
    12.621951219512194,
    10.792682926829267,
    9.512195121951219,
    7.134146341463413,
    5.853658536585366,
])

y = np.array([
    0.7999999999999989,
    1.0434782608695627,
    1.2695652173913032,
    1.530434782608694,
    1.7565217391304326,
    2.0173913043478233,
    2.3652173913043457,
    2.747826086956519,
    3.1652173913043455,
    3.599999999999998,
    4.03478260869565,
    4.452173913043476,
    4.869565217391302,
    5.304347826086953,
    5.739130434782606,
    6.20869565217391,
    6.591304347826084,
    7.060869565217389,
    7.443478260869563,
    7.8608695652173886,
    8.34782608695652,
    8.52173913043478,
])

import numpy as np

xa2 = np.array([
    5.272727272727273,
    7.272727272727273,
    8.363636363636365,
    10.363636363636365,
    11.454545454545455,
    12.181818181818182,
    13.454545454545457,
    14.000000000000002,
    14.181818181818182,
    14.727272727272728,
    14.181818181818182,
    14.181818181818182,
    13.454545454545457,
    12.909090909090908,
    12.363636363636363,
    9.636363636363638,
    8.363636363636365,
])

ya2 = np.array([
    1.046728971962617,
    1.7570093457943923,
    2.205607476635514,
    2.672897196261683,
    3.1214953271028048,
    3.5140186915887854,
    3.981308411214953,
    4.4485981308411215,
    4.841121495327103,
    5.345794392523365,
    5.757009345794392,
    6.261682242990655,
    6.710280373831775,
    7.177570093457945,
    7.588785046728972,
    8.560747663551401,
    9.102803738317757,
])


y_vals = []
u_vals = []

for cell in cells:
    # triângulo tem 3 arestas
    edges = [(cell[0], cell[1]),
             (cell[1], cell[2]),
             (cell[2], cell[0])]

    for n1, n2 in edges:
        x1, y1 = points[n1]
        x2, y2 = points[n2]

        # Checar cruzamento
        if (x1 - x_cut) * (x2 - x_cut) < 0.0:
            alpha = (x_cut - x1) / (x2 - x1)

            y_int = y1 + alpha * (y2 - y1)
            u_int = u_p[n1] + alpha * (u_p[n2] - u_p[n1])

            y_vals.append(y_int)
            u_vals.append(u_int)

# Converter para array
y_vals = np.array(y_vals)
u_vals = np.array(u_vals)

# ======================================================
# 4. Limpeza e ordenação
# ======================================================
# Remover pontos duplicados (arestas compartilhadas)
data = np.column_stack((y_vals, u_vals))
data = np.unique(data, axis=0)

y_vals = data[:, 0]
u_vals = data[:, 1]

order = np.argsort(y_vals)
y_vals = y_vals[order]
u_vals = u_vals[order]

# ======================================================
# 5. Plot
# ======================================================
# print(values_armalys_xs_3d06[::2].shape)
# print(values_armalys_xs_3d06[1::2].shape)

plt.figure(figsize=(4, 10))
plt.scatter(x, y, marker=r'$\circ$', c='orange', label='Armalys et. al.')

# modelo
x_model = u_vals * 30
y_model = y_vals * 10
plt.plot(x_model, y_model, label='Este trabalho')

# barrinhas horizontais
# for xs, xm, y in zip(xa2, x_model, ya2):
#     plt.hlines(y=y, xmin=min(xs, xm), xmax=max(xs, xm),
#                colors='gray', linestyles='dashed', linewidth=1)
plt.xlabel(r"$u$")
plt.ylabel(r"$y$")
# plt.title("Perfil interpolado em x = 1.08")
plt.xlim((-1,30))
plt.ylim((0,10))
plt.grid(True)
plt.tight_layout()
plt.legend()
plt.savefig("comparison_re100_armalys_306.pdf", format='pdf', dpi=300)
