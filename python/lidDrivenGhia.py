import meshio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.tri import LinearTriInterpolator
import sys

Re = float(sys.argv[1])

# Configurações do estilo global do Matplotlib
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "palatino",
    "axes.labelsize": 14,
    "axes.titlesize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12
})

# --- Leitura do arquivo VTK ---
mesh = meshio.read(sys.argv[2])
points = mesh.points
x, y = points[:, 0], points[:, 1]

# --- Dados de velocidade por célula ---
velocity = mesh.cell_data_dict["velocity"]["quad"]
u_cell, v_cell = velocity[:, 0], velocity[:, 1]
vel_mag_cell = np.hypot(u_cell, v_cell)
quads = mesh.cells_dict["quad"]

# --- Conversão: célula → ponto (média nos vértices) ---
num_points = len(points)
u_points = np.zeros(num_points)
v_points = np.zeros(num_points)
vel_mag_points = np.zeros(num_points)
counts = np.zeros(num_points, dtype=int)

for i, quad in enumerate(quads):
    u_points[quad] += u_cell[i]
    v_points[quad] += v_cell[i]
    vel_mag_points[quad] += vel_mag_cell[i]
    counts[quad] += 1

mask = counts > 0
u_points[mask] /= counts[mask]
v_points[mask] /= counts[mask]
vel_mag_points[mask] /= counts[mask]

# --- Criação da triangulação ---
triangles = np.vstack([
    [q[0], q[1], q[2]] for q in quads
] + [
    [q[0], q[2], q[3]] for q in quads
])
triang = tri.Triangulation(x, y, triangles)
triangles = np.array(triangles).reshape(-1, 3)
triang = tri.Triangulation(x, y, triangles)

# --- Interpolação para a grade regular ---
xi = np.linspace(x.min(), x.max(), 100)
yi = np.linspace(y.min(), y.max(), 100)
Xi, Yi = np.meshgrid(xi, yi)

Ui = LinearTriInterpolator(triang, u_points)(Xi, Yi)
Vi = LinearTriInterpolator(triang, v_points)(Xi, Yi)

# --- Carrega dados de Ghia et al. ---
ghia_vertical = np.loadtxt("../docs/ghia_u_vertical_center_line_y.csv", skiprows=1, delimiter=',', dtype=np.float64)
ghia_horizontal = np.loadtxt("../docs/ghia_v_horizontal_center_line_x.csv", skiprows=1, delimiter=',', dtype=np.float64)

# --- Seleciona a coluna onde x = 0.5 ---
col_idx = np.argmin(np.abs(xi - 0.5))
u_center = Ui[:, col_idx]  # perfil vertical de u em x = 0.5
col_idx = np.argmin(np.abs(yi - 0.5))
v_center = Vi[col_idx, :]  # perfil vertical de u em x = 0.5

# --- Plotagem ---
plt.figure(figsize=(8, 8))
plt.plot(u_center, yi, label="Este trabalho", linewidth=2)
if Re == 100:
    plt.scatter(ghia_vertical[:, 1], ghia_vertical[:, 0], label="Ghia et al.", color='black', facecolors='none', marker='o', linewidths=1.5)
elif Re == 400:
    plt.scatter(ghia_vertical[:, 2], ghia_vertical[:, 0], label="Ghia et al.", color='black', facecolors='none', marker='o', linewidths=1.5)
elif Re == 1000:
    plt.scatter(ghia_vertical[:, 3], ghia_vertical[:, 0], label="Ghia et al.", color='black', facecolors='none', marker='o', linewidths=1.5)
plt.xlabel(r"$u$")
plt.ylabel(r"y")
plt.legend()
plt.grid(True, linestyle='--', alpha=0.5)
plt.gca().set_aspect('equal')
plt.tight_layout()
plt.savefig(sys.argv[3], format='pdf', dpi=400)
plt.close()

plt.figure(figsize=(8, 8))
plt.plot(v_center, xi, label="Este trabalho", linewidth=2)
if Re == 100:
    plt.scatter(ghia_horizontal[:, 1], ghia_horizontal[:, 0], label="Ghia et al.", color='black', facecolors='none', marker='o', linewidths=1.5)
elif Re == 400:
    plt.scatter(ghia_horizontal[:, 2], ghia_horizontal[:, 0], label="Ghia et al.", color='black', facecolors='none', marker='o', linewidths=1.5)
elif Re == 1000:
    plt.scatter(ghia_horizontal[:, 3], ghia_horizontal[:, 0], label="Ghia et al.", color='black', facecolors='none', marker='o', linewidths=1.5)
plt.xlabel(r"$v$")
plt.ylabel(r"y")
plt.legend()
plt.grid(True, linestyle='--', alpha=0.5)
plt.gca().set_aspect('equal')
plt.tight_layout()
plt.savefig(sys.argv[4], format='pdf', dpi=400)
plt.close()
