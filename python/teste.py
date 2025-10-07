import pandas as pd
import matplotlib.pyplot as plt

# Leitura do CSV
csv_path = "../outputs/reynolds1000.csv"
df = pd.read_csv(csv_path)

# Valores de Ghia (da imagem)
y_ghia = [
    1.00000, 0.97660, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344,
    0.6172, 0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703,
    0.0625, 0.0547, 0.0000
]

x_ghia = [
    1.00000, 0.65928, 0.57492, 0.51117, 0.46604, 0.33304, 0.18719,
    0.05702, -0.06080, -0.10648, -0.27805, -0.38289, -0.29730, -0.22220,
    -0.20196, -0.18109, 0.00000
]

# Plot principal
plt.figure(figsize=(6,8))
plt.plot(df['velocity_0'], df['Points_1'], label='Simulação')

# Pontos de referência de Ghia
plt.scatter(x_ghia, y_ghia, color='red', label='Ghia', zorder=5)

plt.xlabel('velocity_0')
plt.ylabel('Points_1 (y)')
plt.title('Comparação com Ghia')
plt.legend()
plt.grid(True)
plt.show()
