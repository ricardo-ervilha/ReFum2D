import meshio
import numpy as np
import matplotlib.pyplot as plt

# ============================
# Parâmetros
# ============================
vtk_file = "../outputs/flow_over_cylinder/Flow Over a cylinder_Re=25.vtk"

y_center = 0.205        # linha horizontal passando pelo centro do cilindro
y_tol = 0.002           # tolerância para capturar células próximas à linha
x_start = 0.0
x_end = 10.0
n_bins = 300            # resolução em x

# ============================
# Leitura do VTK
# ============================
mesh = meshio.read(vtk_file)

points = mesh.points[:, :2]
cells = mesh.cells_dict["triangle"]

# componente u da velocidade (por célula)
u = mesh.cell_data["velocity"][0][:, 0]

# ============================
# Centroide das células
# ============================
cell_centers = points[cells].mean(axis=1)
x_c = cell_centers[:, 0]
y_c = cell_centers[:, 1]

# ============================
# Seleção: linha horizontal no centro do cilindro
# ============================
mask = (
    (np.abs(y_c - y_center) <= y_tol) &
    (x_c >= x_start) &
    (x_c <= x_end)
)

x_sel = x_c[mask]
u_sel = u[mask]

# ============================
# Média de u por bins em x
# ============================
bins = np.linspace(x_start, x_end, n_bins)
x_mid = 0.5 * (bins[:-1] + bins[1:])
u_mean = np.full(len(x_mid), np.nan)

for i in range(len(x_mid)):
    in_bin = (x_sel >= bins[i]) & (x_sel < bins[i+1])
    if np.any(in_bin):
        u_mean[i] = np.mean(u_sel[in_bin])

# Remover NaNs
mask_valid = ~np.isnan(u_mean)
x_mid = x_mid[mask_valid]
u_mean = u_mean[mask_valid]

# ============================
# Encontrar reattachment (u < 0 → u > 0)
# ============================
reattachment_x = None

for i in range(len(u_mean) - 1):
    if u_mean[i] < 0 and u_mean[i+1] > 0:
        x1, x2 = x_mid[i], x_mid[i+1]
        u1, u2 = u_mean[i], u_mean[i+1]
        reattachment_x = x1 - u1 * (x2 - x1) / (u2 - u1)
        break

if reattachment_x is None:
    print("❌ Reattachment não encontrado")
else:
    print(f"✅ Reattachment point x = {reattachment_x:.4f}")

# ============================
# Plot
# ============================
plt.figure(figsize=(9,4))
plt.plot(x_mid, u_mean, label=r'$u(x, y=0.205)$')
plt.axhline(0.0, linestyle='--', color='k')

if reattachment_x is not None:
    plt.axvline(reattachment_x, linestyle='--', color='r',
                label='Reattachment')

plt.xlabel("x")
plt.ylabel("u")
plt.title("Reattachment – Flow over a cylinder (linha central)")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
