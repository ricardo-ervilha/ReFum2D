import numpy as np
import matplotlib.pyplot as plt

ncells = np.array([5,10,15,20,25,30,35])
errors = np.array([0.32, 0.086, 0.083, 0.022, 0.015, 0.01, 0.00749786])

A = np.vstack([np.log(ncells), np.ones(len(ncells))]).T
a, b = np.linalg.lstsq(A, np.log(errors), rcond=None)[0]

plt.scatter(ncells, errors)
plt.title(f"m = {a}")
plt.semilogx()
plt.semilogy()
plt.show()