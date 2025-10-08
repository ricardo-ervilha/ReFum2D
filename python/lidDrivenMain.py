import meshio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.tri import LinearTriInterpolator
import sys

Re = sys.argv[1]

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
xi = np.linspace(x.min(), x.max(), 300)
yi = np.linspace(y.min(), y.max(), 300)
Xi, Yi = np.meshgrid(xi, yi)

Ui = LinearTriInterpolator(triang, u_points)(Xi, Yi)
Vi = LinearTriInterpolator(triang, v_points)(Xi, Yi)

# --- Plotagem ---
plt.figure(figsize=(9, 9))

# Campo de magnitude da velocidade
tpc = plt.tricontourf(triang, vel_mag_points, levels=80, cmap='coolwarm', alpha=0.9)

# Linhas de corrente
plt.streamplot(xi, yi, Ui, Vi, color='k', linewidth=0.8, density=2.0, arrowsize=0.6)

# Barra de cores e rótulos
cbar = plt.colorbar(tpc, orientation='horizontal', pad=0.1, shrink=0.5)
cbar.set_label(r"$|\vec{u}|$", fontsize=14)
cbar.ax.tick_params(labelsize=12)

plt.xlabel("x")
plt.ylabel("y")
plt.title(fr"Streamlines — Lid-driven cavity — Re={Re}", fontsize=16, fontweight='bold')
plt.xlim(min(x), max(x)) 
plt.ylim(min(y), max(y)) 
plt.gca().set_aspect('equal') 
plt.tight_layout()
plt.savefig(sys.argv[3], format='pdf', dpi=400)
plt.close()