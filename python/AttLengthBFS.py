import meshio
import numpy as np
import matplotlib.pyplot as plt

# ============================
# Parâmetros
# ============================
vtk_file = "../outputs/backward_facing_step_Re=100.vtk"
x_step = 2.0
y_max_wall = 0.05      # faixa próxima da parede (ajuste se necessário)
n_bins = 200           # resolução em x

# ============================
# Leitura do VTK
# ============================
mesh = meshio.read(vtk_file)

points = mesh.points[:, :2]
cells = mesh.cells_dict["triangle"]

u = mesh.cell_data["velocity"][0][:, 0]

# ============================
# Centroide das células
# ============================
cell_centers = points[cells].mean(axis=1)
x_c = cell_centers[:, 0]
y_c = cell_centers[:, 1]

# ============================
# Seleção: região próxima à parede inferior pós-degrau
# ============================
mask = (
    (y_c >= 0.0) &
    (y_c <= y_max_wall) &
    (x_c > x_step)
)

x_sel = x_c[mask]
u_sel = u[mask]

# ============================
# Média de u por bins em x
# ============================
bins = np.linspace(x_sel.min(), x_sel.max(), n_bins)
x_mid = 0.5 * (bins[:-1] + bins[1:])
u_mean = np.zeros(len(x_mid))

for i in range(len(x_mid)):
    in_bin = (x_sel >= bins[i]) & (x_sel < bins[i+1])
    if np.any(in_bin):
        u_mean[i] = np.mean(u_sel[in_bin])
    else:
        u_mean[i] = np.nan

# Remover NaNs
mask_valid = ~np.isnan(u_mean)
x_mid = x_mid[mask_valid]
u_mean = u_mean[mask_valid]

# ============================
# Encontrar reattachment (cruzamento de zero)
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
    Lr = reattachment_x - x_step
    print(f"✅ Reattachment point x = {reattachment_x:.4f}")
    print(f"✅ Reattachment length Lr = {Lr*2:.4f}")

# ============================
# Plot
# ============================
plt.figure(figsize=(8,4))
plt.plot(x_mid, u_mean, label=r'$\langle u \rangle$ próximo à parede')
plt.axhline(0.0, linestyle='--', color='k')

if reattachment_x:
    plt.axvline(reattachment_x, linestyle='--', label='Reattachment')

plt.xlabel("x")
plt.ylabel("u")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
