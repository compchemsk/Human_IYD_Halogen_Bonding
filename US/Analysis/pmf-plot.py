import numpy as np
import matplotlib.pyplot as plt

# Load data
data1 = np.loadtxt("fe_c1_it2.dat") #Pathway 1 file
data2 = np.loadtxt("fe_c2_it2.dat") #Pathway 2 file
data3 = np.loadtxt("fe_c3_it2.dat") #Pathway 3 file

x1 = data1[:, 0]
y1 = data1[:, 1]
yerr1 = data1[:, 2]

x2 = data2[:, 0]
y2 = data2[:, 1]
yerr2 = data2[:, 2]

x3 = data3[:, 0]
y3 = data3[:, 1]
yerr3 = data3[:, 2]

# Plot with error bars
plt.figure()
plt.errorbar(x1, y1, yerr=yerr1, linewidth=2, color='black', fmt='o-', capsize=3, label='Pathway A')
plt.errorbar(x2, y2, yerr=yerr2, linewidth=2, color='red', fmt='o-', capsize=3, label='Pathway B')
plt.errorbar(x3, y3, yerr=yerr3, linewidth=2, color='blue', fmt='o-', capsize=3, label='Pathway C')

# Labels
plt.xlabel('Progress Coordinate (\u212B)', fontsize=24, fontweight='bold')
plt.ylabel("PMF (kcal/mol)", fontsize=24, fontweight='bold')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(loc='lower right', frameon=False, fontsize=18)
plt.xlim(3.8,19.8)
plt.ylim(0,30)

# Improve layout
plt.tight_layout()
#plt.axvline(x=6.6,linestyle='--',color='grey')

# Axis styling
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_color('black')
    spine.set_linewidth(2)

# Save and show
plt.savefig('pmf-us.png')
plt.show()