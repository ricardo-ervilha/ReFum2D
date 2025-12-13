import meshio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.tri import LinearTriInterpolator
import sys

# ============================================================
#  Entrada
# ============================================================
Re = float(sys.argv[1])
vtk_file = sys.argv[2]

# ============================================================
#  Estilo global Matplotlib
# ============================================================
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "palatino",
    "axes.labelsize": 14,
    "axes.titlesize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12
})

# ============================================================
#  Leitura do VTK
# ============================================================
mesh = meshio.read(vtk_file)

points = mesh.points
x, y = points[:, 0], points[:, 1]

cells = mesh.cells_dict
cell_data = mesh.cell_data_dict["velocity"]

# ============================================================
#  Coleta de velocidades e conectividade
# ============================================================
u_cell = []
v_cell = []
triangles = []

cell_idx = 0
for cell_type, elems in cells.items():

    vel = cell_data[cell_type]

    for i, elem in enumerate(elems):
        u_cell.append(vel[i, 0])
        v_cell.append(vel[i, 1])

        if cell_type == "quad":
            triangles.append([elem[0], elem[1], elem[2]])
            triangles.append([elem[0], elem[2], elem[3]])
        elif cell_type == "triangle":
            triangles.append(elem)

u_cell = np.array(u_cell)
v_cell = np.array(v_cell)

triangles = np.array(triangles)

# ============================================================
#  Conversão célula → ponto (média nos vértices)
# ============================================================
num_points = len(points)
u_points = np.zeros(num_points)
v_points = np.zeros(num_points)
counts = np.zeros(num_points, dtype=int)

cell_idx = 0
for cell_type, elems in cells.items():
    for elem in elems:
        for node in elem:
            u_points[node] += u_cell[cell_idx]
            v_points[node] += v_cell[cell_idx]
            counts[node] += 1
        cell_idx += 1

mask = counts > 0
u_points[mask] /= counts[mask]
v_points[mask] /= counts[mask]

# ============================================================
#  Triangulação
# ============================================================
triang = tri.Triangulation(x, y, triangles)

# ============================================================
#  Interpolação para grade regular
# ============================================================
xi = np.linspace(x.min(), x.max(), 100)
yi = np.linspace(y.min(), y.max(), 100)
Xi, Yi = np.meshgrid(xi, yi)

Ui = LinearTriInterpolator(triang, u_points)(Xi, Yi)
Vi = LinearTriInterpolator(triang, v_points)(Xi, Yi)

# ============================================================
#  Dados de Ghia et al.
# ============================================================
ghia_vertical = np.loadtxt(
    "../docs/ghia_u_vertical_center_line_y.csv",
    skiprows=1, delimiter=',', dtype=np.float64
)

ghia_horizontal = np.loadtxt(
    "../docs/ghia_v_horizontal_center_line_x.csv",
    skiprows=1, delimiter=',', dtype=np.float64
)

# ============================================================
#  Perfis centrais
# ============================================================
col_x = np.argmin(np.abs(xi - 0.5))
u_center = Ui[:, col_x]

col_y = np.argmin(np.abs(yi - 0.5))
v_center = Vi[col_y, :]

# ============================================================
#  Plot — perfil vertical de u
# ============================================================
plt.figure(figsize=(8, 8))
plt.plot(u_center, yi, label="Este trabalho", linewidth=2)

if Re == 100:
    plt.scatter(ghia_vertical[:, 1], ghia_vertical[:, 0],
                label="Ghia et al.", facecolors='none',
                edgecolors='black', marker='o', linewidths=1.5)
elif Re == 400:
    plt.scatter(ghia_vertical[:, 2], ghia_vertical[:, 0],
                label="Ghia et al.", facecolors='none',
                edgecolors='black', marker='o', linewidths=1.5)
elif Re == 1000:
    plt.scatter(ghia_vertical[:, 3], ghia_vertical[:, 0],
                label="Ghia et al.", facecolors='none',
                edgecolors='black', marker='o', linewidths=1.5)

plt.xlabel(r"$u$")
plt.ylabel(r"$y$")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.gca().set_aspect("equal")
plt.tight_layout()
plt.savefig(f"ghia_vertical_{Re}.pdf", dpi=400)
plt.close()

# ============================================================
#  Plot — perfil horizontal de v
# ============================================================
plt.figure(figsize=(8, 8))
plt.plot(v_center, xi, label="Este trabalho", linewidth=2)

if Re == 100:
    plt.scatter(ghia_horizontal[:, 1], ghia_horizontal[:, 0],
                label="Ghia et al.", facecolors='none',
                edgecolors='black', marker='o', linewidths=1.5)
elif Re == 400:
    plt.scatter(ghia_horizontal[:, 2], ghia_horizontal[:, 0],
                label="Ghia et al.", facecolors='none',
                edgecolors='black', marker='o', linewidths=1.5)
elif Re == 1000:
    plt.scatter(ghia_horizontal[:, 3], ghia_horizontal[:, 0],
                label="Ghia et al.", facecolors='none',
                edgecolors='black', marker='o', linewidths=1.5)

plt.xlabel(r"$v$")
plt.ylabel(r"$y$")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.gca().set_aspect("equal")
plt.tight_layout()
plt.savefig(f"ghia_horizontal_{Re}.pdf", dpi=400)
plt.close()
