import numpy as np
import matplotlib.pyplot as plt

N = 30 # num points
# domain limits
x0, xf = 0, 1
y0, yf = x0, xf
h = float((xf - x0)/(N-1))
tol = 1e-8

# parameters
rho = 1.0
u = 1.0
v = 1.0
gamma = 0.1

# face area (2D with unit depth) = h
A_face = h

# face fluxes (convective mass flux through face) -- sign follows velocity direction
F_e = rho * u * A_face
F_w = rho * u * A_face
F_n = rho * v * A_face
F_s = rho * v * A_face

# ---- CORREÇÃO: condutâncias difusivas ----
# D_f = Gamma * A_face / delta = gamma * h / h = gamma
D_e = gamma
D_w = gamma
D_n = gamma
D_s = gamma

# upwind coefficients (standard finite-volume upwind)
a_E = D_e + max(-F_e, 0.0)
a_W = D_w + max( F_w, 0.0)
a_N = D_n + max(-F_n, 0.0)
a_S = D_s + max( F_s, 0.0)

# a_P: soma dos coeficientes dos vizinhos (em escoamento divergência nula isto é suficiente)
a_P = a_E + a_W + a_N + a_S

# domain points (centros dos volumes)
x = np.linspace(x0, xf, N)
y = np.linspace(y0, yf, N)
X, Y = np.meshgrid(x, y, indexing='xy')

# ----- Condições de contorno Dirichlet -----
def phi_w(y):  # x = x0 (esquerda)
    return 1.0
def phi_e(y):  # x = xf (direita)
    return 0.0
def phi_s(x):  # y = y0 (baixo)
    return 0.0
def phi_n(x):  # y = yf (cima)
    return 1.0

# ----- Mapeamento (i,j) -> índice global k -----
def idx(i, j, N):
    # i: linha (y index), j: coluna (x index)
    return i * N + j

# ----- Montagem de A e b -----
nunk = N * N
A = np.zeros((nunk, nunk), dtype=np.float64)
b = np.zeros((nunk,), dtype=np.float64)

for i in range(N):
    for j in range(N):
        p = idx(i, j, N)

        # Nó de contorno? Impõe Dirichlet forte: A_pp = 1, b_p = phi
        if j == 0:
            A[p, p] = 1.0
            b[p] = phi_w(y[i])
            continue
        if j == N - 1:
            A[p, p] = 1.0
            b[p] = phi_e(y[i])
            continue
        if i == 0:
            A[p, p] = 1.0
            b[p] = phi_s(x[j])
            continue
        if i == N - 1:
            A[p, p] = 1.0
            b[p] = phi_n(x[j])
            continue

        # Nó interno
        A[p, p] = a_P

        # vizinhos (todos eles são internos porque este é um nó interno)
        e = idx(i, j + 1, N)
        w = idx(i, j - 1, N)
        n = idx(i + 1, j, N)
        s = idx(i - 1, j, N)

        A[p, e] = -a_E
        A[p, w] = -a_W
        A[p, n] = -a_N
        A[p, s] = -a_S

        # fonte = 0 -> b[p] já é 0


print(A)
print(b)
print("--------------------------------------------------")
# ---- resolve o sistema ----
phi = np.linalg.solve(A, b)

# ---- reorganiza para matriz 2D e plota ----
Phi = phi.reshape((N, N))
plt.figure(figsize=(6,5))
cs = plt.imshow(Phi, origin='lower', cmap='coolwarm')
plt.colorbar(cs)
plt.title("Solução φ (convecção+difusão, upwind)")
plt.xlabel("x")
plt.ylabel("y")
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

# imprime algumas informações
print("Dimensões: A = {}, b = {}".format(A.shape, b.shape))
