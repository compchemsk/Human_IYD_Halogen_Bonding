import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

# Read the whitespace-delimited file
df = pd.read_csv("d-theta-ityr.csv", delim_whitespace=True, header=0) #Contains the data for all pathways

# Normalize trajectory numbers for color mapping
def normalize(series):
    return (series - series.min()) / (series.max() - series.min() + 1e-5)

# Define a function to create a color with varying alpha from a base color
def get_shade(base_color, norm_val):
    rgba = mcolors.to_rgba(base_color)
    darkened = tuple(np.clip(channel * (0.5 + 0.5 * norm_val), 0, 1) for channel in rgba[:3]) + (0.9,)
    return darkened

# Create plot
#plt.figure(figsize=(10, 6))

# Normalize trajectory values for color intensity
norm_trajA = normalize(df['trajA'])
norm_trajB = normalize(df['trajB'])
norm_trajC = normalize(df['trajC'])

# Plot Pathway A (wheat shades)
for i in range(len(df)):
    plt.scatter(df.loc[i, 'distA'], df.loc[i, 'angleA'],
                color=get_shade('wheat', norm_trajA[i]), label='Pathway A' if i == 0 else "", edgecolors='k', linewidths=0.3)

# Plot Pathway B (lightblue shades)
for i in range(len(df)):
    plt.scatter(df.loc[i, 'distB'], df.loc[i, 'angleB'],
                color=get_shade('lightblue', norm_trajB[i]), label='Pathway B' if i == 0 else "", edgecolors='k', linewidths=0.3)

# Plot Pathway C (lightpink shades)
for i in range(len(df)):
    plt.scatter(df.loc[i, 'distC'], df.loc[i, 'angleC'],
                color=get_shade('lightpink', norm_trajC[i]), label='Pathway C' if i == 0 else "", edgecolors='k', linewidths=0.3)

# Labels and legend
# Labels and title
plt.xlabel("d$_{X-A}$ (\u212B)",fontsize=20, fontweight='bold')
plt.ylabel("\u03B8 (degree)",fontsize=20, fontweight='bold')
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
legend = plt.legend(loc='upper right', fontsize=16, frameon=False)
# Color the legend text to match the pathway colors
legend.get_texts()[0].set_color('wheat')      # Pathway A
legend.get_texts()[1].set_color('lightblue')  # Pathway B
legend.get_texts()[2].set_color('lightpink')  # Pathway C
plt.xlim(3.0,4.5)
plt.ylim(140,180)
plt.tight_layout()
plt.show()