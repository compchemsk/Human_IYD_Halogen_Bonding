import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

d1 = np.loadtxt('2d-halo-ener.dat')

# Scatter plot
sc = plt.scatter(
    d1[:,0*10],
    d1[:,1],
    c=d1[:,2],
    cmap="jet",
    edgecolors="k",
    alpha=1,
    vmin=0.0,
    vmax=3.0
)

# Add color bar
cbar = plt.colorbar(sc)
cbar.set_label("Halogen Bond Stabilization Energy \n (kcal/mol)",fontsize=15)

# Labels and title
plt.xlabel("d$_{X-D}$ (\u212B)",fontsize=20)
plt.ylabel("\u03B8 (degree)",fontsize=20)

plt.axhline(y=160, color='red', linestyle='--', linewidth=2)
plt.axvline(x=3.7, color='red', linestyle='--', linewidth=2)

plt.xlim(3.0,4.0)
plt.ylim(140,180)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)


# Show plot
plt.show()
