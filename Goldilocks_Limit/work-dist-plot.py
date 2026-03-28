import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
data = np.loadtxt('work-dist-smd.dat')
colors=['black','red','green','magenta','blue','orange','violet']
labels=['k = 0.5', 'k = 1.0', 'k = 2.0', 'k = 3.0', 'k = 4.0']
plt.figure()
for i in range(2,7):
  sns.kdeplot(data[:, i], color=colors[i % len(colors)])
plt.xlabel('Work (kcal.mol$^{-1}$.\u212B$^{-2}$)', fontsize=20)
plt.ylabel('Counts', fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(labels=labels,frameon=False, fontsize=10)
plt.savefig('work-dist-smd.png')
plt.show()