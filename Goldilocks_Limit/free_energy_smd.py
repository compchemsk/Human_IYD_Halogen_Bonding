import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('work-dist-smd.dat')

list1 = []

for i in range(2, 9):  # Adjust based on your actual data columns
    data1 = np.exp(-data[:, i] / 0.6)
    delG1 = -0.6 * np.log(np.mean(data1))
    list1.append(delG1)

x = np.array([0.1, 0.25, 0.5, 1.0, 2.0, 3.0, 4.0])
list1 = np.array(list1)

plt.plot(x, list1, color='black', marker='*', alpha=0.7)

plt.xlabel('Force constant (kcal.mol$^{-1}$.$\\AA^{-2}$)', fontsize=20)
plt.ylabel('$\\Delta G$ (kcal.mol$^{-1}$)', fontsize=20)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.xlim(0, 4.5)
plt.ylim(80, 200)

plt.savefig('free_energy_smd.png')
plt.show()