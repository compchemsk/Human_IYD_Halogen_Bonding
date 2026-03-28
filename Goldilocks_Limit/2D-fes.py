import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import math 
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
rcParams.update({'font.size': 15})


x, y, z = np.genfromtxt(r'xxx', unpack=True)

plt.figure(1)
tt = plt.tricontourf(x,y,z, cmap="jet")

#ticks = [0.0, -1.0, -2.0, -3.0, -4.0, -5.0]
#cbar = plt.colorbar(tt, ticks=ticks)
#cbar.ax.set_yticklabels([f'{tick}' for tick in ticks])

plt.xlabel('IC1 (d in \u212B)',fontsize=20.0)
plt.ylabel('IC2 (\u03B8)',fontsize=20.0)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.autoscale()
plt.colorbar()
plt.xlim(-40,40)
plt.ylim(0,180)
plt.show()
