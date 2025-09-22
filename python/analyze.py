import matplotlib.pyplot as plt
import numpy as np

Y_values = [10, 20, 30, 40, 50, 100, 200]

errors = []
hs = []

for Y in Y_values:
    suffix = f"{Y}x{Y}"
    filepath = f"./metric_{suffix}.txt"

    try:
        with open(filepath, "r") as f:
            error = float(f.read().strip())
            errors.append(error)
            hs.append((1.0 / Y))  # tamanho da discretização (DIVIDIR POR sqrt(2) quando for triang transfinite)
            print(f"Lido: {filepath} -> erro = {error}")
    except FileNotFoundError:
        print(f"Arquivo {filepath} não encontrado, pulando...")

# Plotando h vs erro
plt.figure(figsize=(8, 6))
plt.loglog(hs, errors, marker="o", linestyle="-", label="Erro numérico")

plt.xlabel("h (tamanho da malha)")
plt.ylabel("Erro")

log_h = np.log(hs)
log_err = np.log(errors)
p, logC = np.polyfit(log_h, log_err, 1)

plt.title(f"Convergência: h vs Erro. Ordem de convergência: {p:.3f}")

plt.grid(True, which="both", linestyle="--", alpha=0.7)
plt.legend()

plt.show()
plt.close()