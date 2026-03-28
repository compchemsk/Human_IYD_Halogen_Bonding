import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Load the data
df = pd.read_csv("hbond_all.dat", delim_whitespace=True, header=None)

# np values (for distance representation)
np_values = [3.8,4.0,4.2,4.4,4.6,4.8,5.0,5.2,5.4,5.6,5.8,6.0,6.2,6.4,6.6,6.8,7.0,7.2,7.4,7.6,7.8,8.0,8.2,8.4,8.6,8.8,9.8,10.8,11.8,12.8,13.8,14.8,15.8,16.8,17.8,18.8]

# Colormap
cmap = cm.get_cmap('tab20c', len(np_values))

# Create a color dictionary for residue colors
residue_color_dict = {}

# 6x6 grid layout
fig, axes = plt.subplots(6, 6, figsize=(36, 28), sharey=True, constrained_layout=True)
axes = axes.flatten()

for i in range(36):
    ax = axes[i]

    x_col = 2 * i
    y_col = 2 * i + 1

    data = df[[x_col, y_col]].dropna()
    x_raw = data[x_col].astype(int)
    y_vals = data[y_col]

    filtered_data = data[y_vals > 5]
    x_raw_filtered = filtered_data[x_col].astype(int)
    y_vals_filtered = filtered_data[y_col]

    if len(y_vals_filtered) == 0:
        print(f"Warning: No data with % H-Bonds > 5 for plot {i+1}, skipping.")
        continue

    x_transformed = []
    x_labels = []
    color_list = []

    for x in x_raw_filtered:
        if x < 221:
            x_new = x + 70
            x_label = str(x_new)
        elif 221 <= x <= 440:
            x_new = x - 150
            x_label = f"{x_new}*"
        elif x == 441:
            x_new = 'FMNhq'
            x_label = 'FMNhq'
        else:
            x_new = x
            x_label = str(x_new)

        if x_new not in residue_color_dict:
            residue_color_dict[x_new] = cmap(len(residue_color_dict) % len(np_values))

        color_list.append(residue_color_dict[x_new])
        x_transformed.append(x_new)
        x_labels.append(x_label)

    bar_positions = list(range(len(x_labels)))

    ax.bar(bar_positions, y_vals_filtered, color=color_list, width=0.4)

    distance_value = np_values[i]
    ax.set_title(f'PC (ζ) = {distance_value} \u212B', fontsize=14)

    ax.set_xticks(bar_positions)
    ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=10)

    ax.tick_params(axis='y', labelsize=12)

    if i % 6 == 0:
        ax.set_ylabel('% H-Bonds', fontsize=14)
    if i >= 30:
        ax.set_xlabel('Residue', fontsize=14)

    ax.set_ylim(0, 100)

plt.show()
