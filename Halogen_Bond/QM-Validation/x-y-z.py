import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('2d-halo-ener.dat')
x = data[:, 0]
y = data[:, 1]
z = data[:, 2]

# Normalize and apply log
z_transformed = np.log1p(z - np.min(z))  # log1p avoids log(0)

# Scale z values for point sizes
z_size_scaled = 20 + 180 * (z - np.min(z)) / (np.max(z) - np.min(z))  # size range: 20 to 200

sc = plt.scatter(
    x,
    y,
    c=z_transformed,             # color by transformed z
    s=z_size_scaled,             # size by original z
    cmap="viridis",
    edgecolors="k",
    alpha=1,
    vmin=0.0,
    vmax=1.0
)

# Add color bar
cbar = plt.colorbar(sc)
cbar.set_label("X-Bonding Strength (kcal/mol)", fontsize=18)

# Labels and title
plt.xlabel("d$_{X-A}$ (\u212B)", fontsize=20)
plt.ylabel("\u03B8 (degree)", fontsize=20)

# Add lines
plt.axhline(y=160, color='red', linestyle='--', linewidth=2)
plt.axvline(x=3.7, color='red', linestyle='--', linewidth=2)

plt.xlim(3.0, 4.0)
plt.ylim(140, 180)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()
