import numpy as np
import matplotlib.pyplot as plt

# --- parâmetros do problema ---
Nx, Ny = 100,100        # número de pontos na malha
x = np.linspace(0, 1, Nx)
y = np.linspace(0, 1, Ny)
dx = x[1] - x[0]
dy = y[1] - y[0]

# Termo fonte
f = 2.0

# Inicializa a matriz de solução
phi = np.zeros((Ny, Nx))

# --- Condições de Dirichlet ---
def phi_dirichlet(x, y):
    return x**2 + y

# --- Condição de Neumann na borda direita ---
def neumann_right(x, y):
    return 2.0 * x  # ∂φ/∂n = 2x

# --- Montagem da matriz do sistema A*phi = b ---
A = np.zeros((Nx*Ny, Nx*Ny))
b = np.zeros(Nx*Ny)

# Índice para converter (i,j) em 1D
def idx(i, j):
    return i + j*Nx

for j in range(Ny):
    for i in range(Nx):
        k = idx(i, j)

        # --- BORDAS ---
        # Esquerda - Dirichlet
        if i == 0:
            A[k, k] = 1.0
            b[k] = phi_dirichlet(x[i], y[j])
            continue
        # Direita - Neumann
        if i == Nx-1:
            A[k, k] = -1.0/dx  # aproximação de ∂φ/∂x
            A[k, k-1] = 1.0/dx
            b[k] = neumann_right(x[i], y[j])
            continue
        # Baixo - Dirichlet
        if j == 0:
            A[k, k] = 1.0
            b[k] = phi_dirichlet(x[i], y[j])
            continue
        # Cima - Dirichlet
        if j == Ny-1:
            A[k, k] = 1.0
            b[k] = phi_dirichlet(x[i], y[j])
            continue

        # --- PONTO INTERNO ---
        A[k, idx(i-1, j)] = 1.0/dx**2
        A[k, idx(i+1, j)] = 1.0/dx**2
        A[k, idx(i, j-1)] = 1.0/dy**2
        A[k, idx(i, j+1)] = 1.0/dy**2
        A[k, k] = -2.0*(1/dx**2 + 1/dy**2)
        b[k] = -f

# --- Resolve o sistema ---
phi_vector = np.linalg.solve(A, b)
phi = phi_vector.reshape((Ny, Nx))

# --- Solução exata ---
X, Y = np.meshgrid(x, y)
phi_exact = X**2 + Y

# --- Plot ---
plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
plt.imshow(phi, origin='lower', cmap='coolwarm', extent=[0,1,0,1])
plt.title("Solução Numérica")
plt.colorbar()

plt.subplot(1,2,2)
plt.imshow(phi_exact, origin='lower', cmap='turbo', extent=[0,1,0,1])
plt.title("Solução Exata")
plt.colorbar()
plt.show()

# --- Erro ---
error = np.max(np.abs(phi - phi_exact))
print("Erro máximo:", error)
