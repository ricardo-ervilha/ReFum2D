import numpy as np
import matplotlib.pyplot as plt

ncells = np.array([5,10,15,20,25,30,35])
errors = np.array([0.11,0.03, 0.0237, 0.008, 0.0052, 0.0036, 0.0026])

A = np.vstack([np.log(ncells), np.ones(len(ncells))]).T
a, b = np.linalg.lstsq(A, np.log(errors), rcond=None)[0]

plt.scatter(ncells, errors)
plt.title(f"m = {a}")
plt.semilogx()
plt.semilogy()
plt.show()