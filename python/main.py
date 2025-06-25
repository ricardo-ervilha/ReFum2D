import numpy as np
import matplotlib.pyplot as plt

ncells = np.array([5,10,15,20,25,30,35])
errors = np.array([0.19,0.0486,0.0358,0.0122,0.00782,0.00543,0.003994])

A = np.vstack([np.log(ncells), np.ones(len(ncells))]).T
a, b = np.linalg.lstsq(A, np.log(errors), rcond=None)[0]

plt.scatter(ncells, errors)
plt.title(f"m = {a}")
plt.semilogx()
plt.semilogy()
plt.show()