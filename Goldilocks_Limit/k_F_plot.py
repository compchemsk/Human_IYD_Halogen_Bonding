import numpy as np
import matplotlib.pyplot as plt
import os

# Your k directories
k_values = [0.1, 0.25, 0.5, 1.0, 2.0, 3.0, 4.0]

plt.figure()

for k in k_values:
    folder = f'k{k}'
    file_path = os.path.join(folder, 'data_force.dat')
    
    data = np.loadtxt(file_path)
    
    x = data[:, 0]   # progress coordinate (distance)
    y = data[:, 1]   # force
    
    plt.plot(x, y, label=f'k = {k}')

plt.xlabel('Reaction Coordinate (Å)', fontsize=20)
plt.ylabel('Force (kcal/mol·Å)', fontsize=20)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(frameon=False, fontsize=10)

plt.savefig('force_vs_progress_coordinate.png')
plt.show()