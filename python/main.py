import meshio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

# * Configurações da matplot
plt.rcParams.update({
    "text.usetex": True,          
    "font.family": "serif",
    "font.serif": ["Computer Modern"],
    "axes.labelsize": 14,
    "axes.titlesize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12
})

# ====================== CÓDIGO GERADO ========================================
mesh = meshio.read("../outputs/benchmarkConvection.vtk")

# Pega os pontos
points = mesh.points
x = points[:, 0]
y = points[:, 1]

# Pega o campo "Temperatura" em CELL_DATA
T_cell = mesh.cell_data_dict["Temperatura"]["quad"]  # quad = tipo da célula

# Cria uma lista de triângulos a partir dos quadriláteros
quads = mesh.cells_dict["quad"]  # cada linha: indices dos 4 pontos
triangles = []
for q in quads:
    # dividir cada quad em 2 triângulos: (0,1,2) e (0,2,3)
    triangles.append([q[0], q[1], q[2]])
    triangles.append([q[0], q[2], q[3]])

triangles = np.array(triangles)

# Duplicar valores de T para os dois triângulos de cada quad
T_tri = np.repeat(T_cell, 2)

# Cria triangulação
triang = tri.Triangulation(x, y, triangles)


# ====================== CÓDIGO GERADO ========================================

plt.figure(figsize=(6,5))
plt.tripcolor(triang, T_tri, shading='flat', cmap='coolwarm')
plt.colorbar(label="u(x,y)")
plt.xlabel("x")
plt.ylabel("y")
plt.tight_layout()
plt.show()
