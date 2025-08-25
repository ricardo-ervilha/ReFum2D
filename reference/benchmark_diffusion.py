import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

N = 100 # num points
# domain limits
x0, xf = 0, 1
y0, yf = x0, xf
gamma = 1.0 # diffusion cte
tol = 1e-8 # gauss-seidel tolerance
h = (xf-x0)/(N-1) # discretization
Af = (xf-x0)*(yf-y0) # area of control volume

# domain points
x = np.linspace(x0, xf, N)
y = np.linspace(y0, yf, N)
X,Y = np.meshgrid(x,y)

def exact(x,y):
    return 100 * x * (1-x) * y * (1-y)

def source(x,y):
    return -200*x*(1-x) - 200*y*(1-y)

def poisson_solver():
    # all bc's are zero. 
    u_old = np.zeros((N,N))
    u_new = np.zeros((N,N))

    error = 1
    iter = 0
    while error > tol:
        # solution with for
        # for i in range(1,N-1):
        #     for j in range(1,N-1):
        #         u_new[i,j] = (u_old[i+1,j] + u_old[i-1,j] + u_old[i,j+1] + u_old[i,j-1] - \
        #                        h**2 * source(X[i,j], Y[i,j]))/4.0
        # solution with vectorization
        u_new[1:-1, 1:-1] = (
            u_old[2:, 1:-1] +   # i+1
            u_old[:-2, 1:-1] +  # i-1
            u_old[1:-1, 2:] +   # j+1
            u_old[1:-1, :-2] -  # j-1
            h**2 * source(X[1:-1, 1:-1], Y[1:-1, 1:-1])
        ) / 4.0
        error = np.linalg.norm(u_new - u_old, ord=np.inf)
        u_new, u_old = u_old, u_new
        iter += 1
        if iter % 100 == 0:
            print(f"Error = {error}")
    
    return u_new

# calls to solve
u = poisson_solver()

# plot solution -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
plt.imshow(u, extent=[0,1,0,1], origin="lower", cmap="coolwarm")

plt.xlabel("x", fontsize=14)
plt.xticks(fontsize=14)

plt.ylabel("y", fontsize=14)
plt.yticks(fontsize=14)

cbar = plt.colorbar()
cbar.set_label("u", fontsize=14)

plt.show()
plt.close()
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=