import numpy as np
import matplotlib.pyplot as plt
import os

# Force constants (directory names)
k_values = [0.1, 0.25, 0.5, 1.0, 2.0, 3.0, 4.0]

plt.figure(figsize=(6,5))

for k in k_values:
    folder = f'k{k}'
    file_path = os.path.join(folder, 'data_PMF.dat')
    
    data = np.loadtxt(file_path)
    
    x = data[:, 0]   # reaction coordinate (distance)
    pmf = data[:, 1] # PMF (free energy)

    # Optional: shift PMF so it starts from zero
    pmf = pmf - pmf[0]

    plt.plot(x, pmf, label=f'k = {k}')

plt.xlabel('Reaction Coordinate (Å)', fontsize=20)
plt.ylabel('PMF (kcal/mol)', fontsize=20)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(frameon=False, fontsize=10)
plt.tight_layout()
plt.savefig('pmf_vs_progress_coordinate.png', dpi=300)
plt.show()