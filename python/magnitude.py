import meshio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.tri import LinearTriInterpolator

# Configurações da Matplotlib
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "palatino",
    "font.serif": ["Computer Modern"],
    # tamanhos
    "axes.labelsize": 14,
    "axes.titlesize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12
})

# --- Leitura da Malha ---
mesh = meshio.read("../outputs/v_vector.vtk")

# Coordenadas dos pontos da malha
points = mesh.points
x = points[:, 0]
y = points[:, 1]
num_points = len(points)

# --- Extração dos Dados das Células ---
# Pega os dados de velocidade associados às células quadrilaterais
velocity_cell = mesh.cell_data_dict["velocity"]["quad"]
u_cell = velocity_cell[:, 0]
v_cell = velocity_cell[:, 1]
vel_mag_cell = np.sqrt(u_cell**2 + v_cell**2)

# Pega a conectividade dos quadriláteros
quads = mesh.cells_dict["quad"]

# --- Conversão de Dados de Célula para Ponto (Nodal) ---
# Inicializa arrays para armazenar os valores nos pontos
u_points = np.zeros(num_points)
v_points = np.zeros(num_points)
vel_mag_points = np.zeros(num_points)
point_counts = np.zeros(num_points, dtype=int) # Contador para a média

# Itera sobre cada célula (quad) para distribuir seus valores para seus pontos
for i, quad_cell in enumerate(quads):
    for point_index in quad_cell:
        u_points[point_index] += u_cell[i]
        v_points[point_index] += v_cell[i]
        vel_mag_points[point_index] += vel_mag_cell[i]
        point_counts[point_index] += 1

# Calcula a média para cada ponto dividindo pela contagem
# Evita divisão por zero para pontos que não pertencem a nenhuma célula
non_zero_counts = np.where(point_counts > 0)
u_points[non_zero_counts] /= point_counts[non_zero_counts]
v_points[non_zero_counts] /= point_counts[non_zero_counts]
vel_mag_points[non_zero_counts] /= point_counts[non_zero_counts]

# --- Criação da Triangulação para Plotagem ---
# Converte os quads em triângulos para a visualização
triangles = []
for q in quads:
    triangles.append([q[0], q[1], q[2]])
    triangles.append([q[0], q[2], q[3]])
triangles = np.array(triangles)
triang = tri.Triangulation(x, y, triangles)

# --- Geração do Gráfico ---
plt.figure(figsize=(9, 9))

# 1. Contorno preenchido suave da magnitude (agora usando dados nodais)
tpc = plt.tricontourf(triang, vel_mag_points, levels=80, cmap='coolwarm', alpha=0.9)

# 2. Streamlines (linhas de corrente)
# Cria uma grade regular para interpolar os dados de velocidade
xi = np.linspace(min(x), max(x), 300)
yi = np.linspace(min(y), max(y), 300)
Xi, Yi = np.meshgrid(xi, yi)

# Interpola os componentes u e v (nodais) para a grade regular
interp_u = tri.LinearTriInterpolator(triang, u_points)
Ui = interp_u(Xi, Yi)
interp_v = tri.LinearTriInterpolator(triang, v_points)
Vi = interp_v(Xi, Yi)

# Plota as streamlines
plt.streamplot(
    xi, yi, Ui, Vi,
    color='k', linewidth=0.8, density=2.0, arrowsize=0.6
)

# --- Configurações Finais do Gráfico ---
cbar = plt.colorbar(tpc, orientation='horizontal', pad=0.1, shrink=0.5)
cbar.set_label(r"$\mid$u$\mid$", fontsize=14)
cbar.ax.tick_params(labelsize=12)

plt.xlabel("x", fontsize=14)
plt.ylabel("y", fontsize=14)
plt.title(r"Streamlines — Lid-driven cavity — Re=1000", fontsize=16, fontweight='bold')

plt.xlim(min(x), max(x))
plt.ylim(min(y), max(y))
plt.gca().set_aspect('equal')

plt.tight_layout()
plt.show()