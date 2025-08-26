import numpy as np
import matplotlib.pyplot as plt
from numba import njit, prange

plt.rcParams['text.usetex'] = True

N = 30 # num points
# domain limits
x0, xf = 0, 1
y0, yf = x0, xf
h = np.float64((xf - x0)/(N-1))
tol = 1e-8

# parameters
rho = 1
u = 1.0
v = 1.0
gamma = 1e-2

F_e = rho * u * h
F_w = rho * u * h
F_n = rho * v * h
F_s = rho * v * h

D_e = gamma / h
D_w = gamma / h
D_n = gamma / h
D_s = gamma / h

a_E = D_e + max(0, -F_e)
a_W = D_w + max(F_w, 0)
a_N = D_n + max(0, -F_n)
a_S = D_s + max(F_s, 0)
a_P = a_E + a_W + a_N + a_S

print(F_e)
print(F_n)
print(D_e)
print(a_E)
print(a_W)
print(a_S)
print(a_P)


print(f"Peclet number = {rho*u/(gamma/h)}")

# domain points
x = np.linspace(x0, xf, N)
y = np.linspace(y0, yf, N)
X,Y = np.meshgrid(x,y)

@njit(parallel=True)
def cd_solver():
    iter = 0
    error = 1.0

    u_old = np.zeros((N,N))
    # impose BCs
    u_old[0,:] = 1.0   # top
    u_old[:,0] = 1.0   # left

    u_new = u_old.copy()

    
    # print(f"Peclet number = {np.sqrt(u**2 + v**2) * (xf - x0) / gamma}")

    while error > tol:
        for i in prange(1, N-1):
            for j in range(1, N-1):
                u_new[i,j] = (a_W*u_old[i, j-1] +
                              a_E*u_old[i, j+1] +
                              a_N*u_old[i-1, j] +
                              a_S*u_old[i+1, j]) / a_P

        error = np.linalg.norm(u_new - u_old, ord=np.inf)
        u_old[:,:] = u_new 
        iter += 1

    return u_new


# calls to solve
u = cd_solver()
print(u)

# plot solution -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
plt.imshow(u, extent=[0,1,0,1], cmap="coolwarm")

plt.xlabel("x", fontsize=14)
plt.xticks(fontsize=14)

plt.ylabel("y", fontsize=14)
plt.yticks(fontsize=14)

cbar = plt.colorbar()
cbar.set_label("u", fontsize=14)

plt.show()
plt.close()
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=